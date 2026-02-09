"""Local desktop app for running and visualizing Xenium analysis outputs."""

from __future__ import annotations

import json
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from PySide6 import QtCore, QtGui, QtWidgets

try:
    from PySide6.QtWebEngineWidgets import QWebEngineView
    WEB_AVAILABLE = True
except Exception:
    QWebEngineView = None
    WEB_AVAILABLE = False


ROOT_DIR = Path(__file__).resolve().parents[1]
RECENT_PATH = Path.home() / ".spatial-analysis-for-dummies" / "recent.json"


@dataclass
class RecentProject:
    data_dir: str
    out_dir: str
    karospace_html: Optional[str]
    last_used: str

    def label(self) -> str:
        return f"{self.data_dir} -> {self.out_dir}"


def _load_recent() -> List[RecentProject]:
    if not RECENT_PATH.exists():
        return []
    try:
        payload = json.loads(RECENT_PATH.read_text())
    except json.JSONDecodeError:
        return []
    projects: List[RecentProject] = []
    for item in payload.get("projects", []):
        projects.append(
            RecentProject(
                data_dir=item.get("data_dir", ""),
                out_dir=item.get("out_dir", ""),
                karospace_html=item.get("karospace_html"),
                last_used=item.get("last_used", ""),
            )
        )
    return projects


def _save_recent(projects: List[RecentProject]) -> None:
    RECENT_PATH.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "projects": [
            {
                "data_dir": p.data_dir,
                "out_dir": p.out_dir,
                "karospace_html": p.karospace_html,
                "last_used": p.last_used,
            }
            for p in projects
        ]
    }
    RECENT_PATH.write_text(json.dumps(payload, indent=2))


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("Spatial Analysis for Dummies")
        self.resize(1200, 800)

        self.process: Optional[QtCore.QProcess] = None
        self._plot_processes: List[QtCore.QProcess] = []
        self.current_out_dir: Optional[Path] = None
        self.current_karospace_html: Optional[Path] = None
        self.recent_projects: List[RecentProject] = _load_recent()

        self._build_ui()
        self._populate_recent()

    def _build_ui(self) -> None:
        root = QtWidgets.QWidget()
        layout = QtWidgets.QHBoxLayout(root)
        layout.setContentsMargins(8, 8, 8, 8)

        self.recent_list = QtWidgets.QListWidget()
        self.recent_list.setMaximumWidth(320)
        self.recent_list.itemSelectionChanged.connect(self._on_recent_selected)

        recent_box = QtWidgets.QGroupBox("Recent Projects")
        recent_layout = QtWidgets.QVBoxLayout(recent_box)
        recent_layout.addWidget(self.recent_list)

        self.load_recent_btn = QtWidgets.QPushButton("Load Outputs")
        self.load_recent_btn.clicked.connect(self._load_selected_recent)
        recent_layout.addWidget(self.load_recent_btn)

        layout.addWidget(recent_box)

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self._build_run_tab(), "Run")
        self.tabs.addTab(self._build_qc_tab(), "QC")
        self.tabs.addTab(self._build_spatial_tab(), "Spatial")
        self.tabs.addTab(self._build_umap_tab(), "UMAP")
        self.tabs.addTab(self._build_compartment_tab(), "Compartments")
        layout.addWidget(self.tabs)

        self.setCentralWidget(root)

    def _build_run_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)

        form = QtWidgets.QFormLayout()

        self.data_dir_edit = QtWidgets.QLineEdit()
        self.data_dir_edit.setPlaceholderText("/path/to/dataset")
        data_btn = QtWidgets.QPushButton("Browse")
        data_btn.clicked.connect(lambda: self._choose_dir(self.data_dir_edit))
        data_row = QtWidgets.QHBoxLayout()
        data_row.addWidget(self.data_dir_edit)
        data_row.addWidget(data_btn)
        form.addRow("Data dir", data_row)

        self.out_dir_edit = QtWidgets.QLineEdit()
        self.out_dir_edit.setPlaceholderText("/path/to/output")
        out_btn = QtWidgets.QPushButton("Browse")
        out_btn.clicked.connect(lambda: self._choose_dir(self.out_dir_edit))
        out_row = QtWidgets.QHBoxLayout()
        out_row.addWidget(self.out_dir_edit)
        out_row.addWidget(out_btn)
        form.addRow("Output dir", out_row)

        self.run_prefix_edit = QtWidgets.QLineEdit("output-")
        form.addRow("Run prefix", self.run_prefix_edit)

        layout.addLayout(form)

        options_box = QtWidgets.QGroupBox("Optional Steps")
        options_layout = QtWidgets.QVBoxLayout(options_box)

        self.mana_check = QtWidgets.QCheckBox("Enable MANA weighted aggregation")
        self.mana_layers = QtWidgets.QSpinBox()
        self.mana_layers.setRange(1, 10)
        self.mana_layers.setValue(3)
        self.mana_hop_decay = QtWidgets.QDoubleSpinBox()
        self.mana_hop_decay.setRange(0.0, 1.0)
        self.mana_hop_decay.setSingleStep(0.05)
        self.mana_hop_decay.setValue(0.5)
        self.mana_kernel = QtWidgets.QComboBox()
        self.mana_kernel.addItems(["exponential", "inverse", "gaussian", "none"])

        mana_form = QtWidgets.QFormLayout()
        mana_form.addRow(self.mana_check)
        mana_form.addRow("Layers", self.mana_layers)
        mana_form.addRow("Hop decay", self.mana_hop_decay)
        mana_form.addRow("Distance kernel", self.mana_kernel)
        options_layout.addLayout(mana_form)

        self.karospace_check = QtWidgets.QCheckBox("Export KaroSpace HTML")
        self.karospace_path_edit = QtWidgets.QLineEdit()
        self.karospace_path_edit.setPlaceholderText("/path/to/karospace.html")
        karospace_btn = QtWidgets.QPushButton("Browse")
        karospace_btn.clicked.connect(self._choose_karospace_path)
        karospace_row = QtWidgets.QHBoxLayout()
        karospace_row.addWidget(self.karospace_path_edit)
        karospace_row.addWidget(karospace_btn)

        options_layout.addWidget(self.karospace_check)
        options_layout.addLayout(karospace_row)

        layout.addWidget(options_box)

        button_row = QtWidgets.QHBoxLayout()
        self.run_btn = QtWidgets.QPushButton("Run Pipeline")
        self.run_btn.clicked.connect(self._run_pipeline)
        self.load_outputs_btn = QtWidgets.QPushButton("Load Outputs Only")
        self.load_outputs_btn.clicked.connect(self._load_outputs_only)
        button_row.addWidget(self.run_btn)
        button_row.addWidget(self.load_outputs_btn)

        layout.addLayout(button_row)

        self.log_view = QtWidgets.QPlainTextEdit()
        self.log_view.setReadOnly(True)
        layout.addWidget(self.log_view, stretch=1)

        return widget

    def _build_qc_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)

        self.qc_scroll = QtWidgets.QScrollArea()
        self.qc_scroll.setWidgetResizable(True)
        self.qc_container = QtWidgets.QWidget()
        self.qc_layout = QtWidgets.QVBoxLayout(self.qc_container)
        self.qc_layout.addStretch(1)
        self.qc_scroll.setWidget(self.qc_container)

        layout.addWidget(self.qc_scroll)
        return widget

    def _build_spatial_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)

        if WEB_AVAILABLE:
            self.spatial_view = QWebEngineView()
            layout.addWidget(self.spatial_view)
        else:
            self.spatial_view = None
            self.spatial_fallback_label = QtWidgets.QLabel(
                "Qt WebEngine not available. Spatial viewer will open in your browser."
            )
            self.spatial_fallback_label.setAlignment(QtCore.Qt.AlignCenter)
            layout.addWidget(self.spatial_fallback_label)
            self.spatial_open_btn = QtWidgets.QPushButton("Open KaroSpace in Browser")
            self.spatial_open_btn.clicked.connect(self._open_karospace_external)
            layout.addWidget(self.spatial_open_btn)

        return widget

    def _build_umap_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)

        self.umap_label = QtWidgets.QLabel("No UMAP image loaded.")
        self.umap_label.setAlignment(QtCore.Qt.AlignCenter)
        self.umap_label.setMinimumHeight(400)
        layout.addWidget(self.umap_label)

        self.umap_generate_btn = QtWidgets.QPushButton("Generate UMAP Plot")
        self.umap_generate_btn.clicked.connect(self._generate_umap_plot)
        layout.addWidget(self.umap_generate_btn)
        return widget

    def _build_compartment_tab(self) -> QtWidgets.QWidget:
        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout(widget)

        self.compartment_label = QtWidgets.QLabel("No compartment map loaded.")
        self.compartment_label.setAlignment(QtCore.Qt.AlignCenter)
        self.compartment_label.setMinimumHeight(400)
        layout.addWidget(self.compartment_label)

        self.compartment_generate_btn = QtWidgets.QPushButton("Generate Compartment Map")
        self.compartment_generate_btn.clicked.connect(self._generate_compartment_map)
        layout.addWidget(self.compartment_generate_btn)
        return widget

    def _choose_dir(self, line_edit: QtWidgets.QLineEdit) -> None:
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Select directory")
        if path:
            line_edit.setText(path)

    def _choose_karospace_path(self) -> None:
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Select KaroSpace HTML",
            "karospace.html",
            "HTML files (*.html)"
        )
        if path:
            self.karospace_path_edit.setText(path)

    def _log(self, message: str) -> None:
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_view.appendPlainText(f"[{timestamp}] {message}")

    def _run_pipeline(self) -> None:
        data_dir = self.data_dir_edit.text().strip()
        out_dir = self.out_dir_edit.text().strip()
        if not data_dir or not out_dir:
            QtWidgets.QMessageBox.warning(self, "Missing paths", "Please set data and output directories.")
            return

        args = [
            sys.executable,
            str(ROOT_DIR / "run_xenium_analysis.py"),
            "--data-dir",
            data_dir,
            "--out-dir",
            out_dir,
            "--run-prefix",
            self.run_prefix_edit.text().strip() or "output-",
        ]

        if self.mana_check.isChecked():
            args += [
                "--mana-aggregate",
                "--mana-n-layers",
                str(self.mana_layers.value()),
                "--mana-hop-decay",
                str(self.mana_hop_decay.value()),
                "--mana-distance-kernel",
                self.mana_kernel.currentText(),
            ]

        if self.karospace_check.isChecked():
            karospace_path = self.karospace_path_edit.text().strip()
            if not karospace_path:
                karospace_path = str(Path(out_dir) / "karospace.html")
                self.karospace_path_edit.setText(karospace_path)
            args += ["--karospace-html", karospace_path]

        self._log("Starting pipeline...")
        self.run_btn.setEnabled(False)
        self.process = QtCore.QProcess(self)
        self.process.setProgram(args[0])
        self.process.setArguments(args[1:])
        self.process.setWorkingDirectory(str(ROOT_DIR))
        self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._on_process_output)
        self.process.finished.connect(self._on_process_finished)
        self.process.start()

    def _on_process_output(self) -> None:
        if not self.process:
            return
        text = self.process.readAllStandardOutput().data().decode("utf-8", errors="ignore")
        for line in text.splitlines():
            if line.strip():
                self._log(line)

    def _on_process_finished(self, exit_code: int, _status: QtCore.QProcess.ExitStatus) -> None:
        self.run_btn.setEnabled(True)
        self._log(f"Pipeline finished (exit code {exit_code}).")
        out_dir = self.out_dir_edit.text().strip()
        if out_dir:
            self._load_outputs(Path(out_dir))
            self._update_recent()

    def _load_outputs_only(self) -> None:
        out_dir = self.out_dir_edit.text().strip()
        if not out_dir:
            out_dir = QtWidgets.QFileDialog.getExistingDirectory(self, "Select output directory")
            if not out_dir:
                return
            self.out_dir_edit.setText(out_dir)
        self._load_outputs(Path(out_dir))
        self._update_recent()

    def _load_outputs(self, out_dir: Path) -> None:
        self.current_out_dir = out_dir
        self._load_qc_images(out_dir)
        self._load_karospace(out_dir)
        self._load_umap_image(out_dir)
        self._load_compartment_image(out_dir)

    def _load_qc_images(self, out_dir: Path) -> None:
        for i in reversed(range(self.qc_layout.count())):
            item = self.qc_layout.takeAt(i)
            if item and item.widget():
                item.widget().deleteLater()
        qc_dir = out_dir / "xenium_qc"
        if not qc_dir.exists():
            self.qc_layout.addWidget(QtWidgets.QLabel("No QC outputs found."))
            self.qc_layout.addStretch(1)
            return

        images = sorted(qc_dir.glob("*.png"))
        if not images:
            self.qc_layout.addWidget(QtWidgets.QLabel("No QC images found."))
            self.qc_layout.addStretch(1)
            return

        for img_path in images:
            label = QtWidgets.QLabel()
            label.setAlignment(QtCore.Qt.AlignCenter)
            pixmap = QtGui.QPixmap(str(img_path))
            if not pixmap.isNull():
                label.setPixmap(pixmap.scaledToWidth(900, QtCore.Qt.SmoothTransformation))
            else:
                label.setText(str(img_path))
            self.qc_layout.addWidget(label)

        self.qc_layout.addStretch(1)

    def _load_karospace(self, out_dir: Path) -> None:
        karospace_path = None
        candidate = out_dir / "karospace.html"
        if candidate.exists():
            karospace_path = candidate
        elif self.karospace_path_edit.text().strip():
            path = Path(self.karospace_path_edit.text().strip())
            if path.exists():
                karospace_path = path

        self.current_karospace_html = karospace_path
        if WEB_AVAILABLE and self.spatial_view is not None:
            if karospace_path:
                self.spatial_view.load(QtCore.QUrl.fromLocalFile(str(karospace_path)))
            else:
                self.spatial_view.setHtml("<p>No KaroSpace HTML found. Run with export enabled.</p>")
        elif not WEB_AVAILABLE:
            if karospace_path:
                self.spatial_fallback_label.setText(f"KaroSpace HTML: {karospace_path}")
            else:
                self.spatial_fallback_label.setText("No KaroSpace HTML found. Run with export enabled.")

    def _open_karospace_external(self) -> None:
        if not self.current_karospace_html or not self.current_karospace_html.exists():
            QtWidgets.QMessageBox.information(
                self,
                "No KaroSpace HTML",
                "No KaroSpace HTML found. Run with export enabled first.",
            )
            return
        QtGui.QDesktopServices.openUrl(
            QtCore.QUrl.fromLocalFile(str(self.current_karospace_html))
        )

    def _load_umap_image(self, out_dir: Path) -> None:
        plot_path = out_dir / "plots" / "umap.png"
        if plot_path.exists():
            pixmap = QtGui.QPixmap(str(plot_path))
            self.umap_label.setPixmap(pixmap.scaledToWidth(900, QtCore.Qt.SmoothTransformation))
        else:
            self.umap_label.setText("No UMAP image found. Click Generate UMAP Plot.")

    def _load_compartment_image(self, out_dir: Path) -> None:
        plot_path = out_dir / "plots" / "compartments.png"
        if plot_path.exists():
            pixmap = QtGui.QPixmap(str(plot_path))
            self.compartment_label.setPixmap(pixmap.scaledToWidth(900, QtCore.Qt.SmoothTransformation))
        else:
            self.compartment_label.setText("No compartment map found. Click Generate Compartment Map.")

    def _generate_umap_plot(self) -> None:
        if not self.current_out_dir:
            QtWidgets.QMessageBox.warning(self, "Missing output", "Load outputs first.")
            return

        h5ad_path = self.current_out_dir / "data" / "clustered.h5ad"
        if not h5ad_path.exists():
            QtWidgets.QMessageBox.warning(self, "Missing file", "clustered.h5ad not found.")
            return

        output_dir = self.current_out_dir / "plots"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "umap.png"

        args = [
            sys.executable,
            "-m",
            "utils.app_visuals",
            "umap",
            "--h5ad",
            str(h5ad_path),
            "--output",
            str(output_path),
        ]

        self._run_visual_process(args, output_path, self.umap_label)

    def _generate_compartment_map(self) -> None:
        if not self.current_out_dir:
            QtWidgets.QMessageBox.warning(self, "Missing output", "Load outputs first.")
            return

        h5ad_path = self.current_out_dir / "data" / "clustered.h5ad"
        if not h5ad_path.exists():
            QtWidgets.QMessageBox.warning(self, "Missing file", "clustered.h5ad not found.")
            return

        output_dir = self.current_out_dir / "plots"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "compartments.png"

        args = [
            sys.executable,
            "-m",
            "utils.app_visuals",
            "compartments",
            "--h5ad",
            str(h5ad_path),
            "--output",
            str(output_path),
        ]

        self._run_visual_process(args, output_path, self.compartment_label)

    def _run_visual_process(
        self,
        args: List[str],
        output_path: Path,
        target_label: QtWidgets.QLabel,
    ) -> None:
        self._log(f"Generating plot: {output_path.name}")
        proc = QtCore.QProcess(self)
        self._plot_processes.append(proc)
        proc.setProgram(args[0])
        proc.setArguments(args[1:])
        proc.setWorkingDirectory(str(ROOT_DIR))
        proc.setProcessChannelMode(QtCore.QProcess.MergedChannels)

        def _on_finished(exit_code: int, _status: QtCore.QProcess.ExitStatus) -> None:
            if proc in self._plot_processes:
                self._plot_processes.remove(proc)
            if exit_code != 0:
                self._log(f"Plot generation failed: {output_path.name}")
                return
            if output_path.exists():
                pixmap = QtGui.QPixmap(str(output_path))
                target_label.setPixmap(pixmap.scaledToWidth(900, QtCore.Qt.SmoothTransformation))
            self._log(f"Plot ready: {output_path.name}")

        proc.finished.connect(_on_finished)
        proc.start()

    def _populate_recent(self) -> None:
        self.recent_list.clear()
        for project in sorted(self.recent_projects, key=lambda p: p.last_used, reverse=True):
            item = QtWidgets.QListWidgetItem(project.label())
            item.setData(QtCore.Qt.UserRole, project)
            self.recent_list.addItem(item)

    def _on_recent_selected(self) -> None:
        items = self.recent_list.selectedItems()
        if not items:
            return
        project: RecentProject = items[0].data(QtCore.Qt.UserRole)
        self.data_dir_edit.setText(project.data_dir)
        self.out_dir_edit.setText(project.out_dir)
        if project.karospace_html:
            self.karospace_path_edit.setText(project.karospace_html)

    def _load_selected_recent(self) -> None:
        items = self.recent_list.selectedItems()
        if not items:
            return
        project: RecentProject = items[0].data(QtCore.Qt.UserRole)
        self.data_dir_edit.setText(project.data_dir)
        self.out_dir_edit.setText(project.out_dir)
        if project.karospace_html:
            self.karospace_path_edit.setText(project.karospace_html)
        self._load_outputs(Path(project.out_dir))

    def _update_recent(self) -> None:
        data_dir = self.data_dir_edit.text().strip()
        out_dir = self.out_dir_edit.text().strip()
        if not data_dir or not out_dir:
            return
        karospace_html = self.karospace_path_edit.text().strip() or None
        now = datetime.now().isoformat(timespec="seconds")

        filtered = [
            p for p in self.recent_projects
            if not (p.data_dir == data_dir and p.out_dir == out_dir)
        ]
        filtered.append(RecentProject(data_dir, out_dir, karospace_html, now))
        self.recent_projects = filtered
        _save_recent(self.recent_projects)
        self._populate_recent()


def main() -> None:
    app = QtWidgets.QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
