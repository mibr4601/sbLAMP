import sys
import os
import json
import re
import tempfile
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QTabWidget,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QLabel,
    QPushButton,
    QTextEdit,
    QFileDialog,
    QGroupBox,
    QDoubleSpinBox,
    QSpinBox,
    QAbstractSpinBox,
    QScrollArea,
    QCheckBox,
    QFrame,
)
from PyQt5.QtCore import Qt, QProcess, pyqtSignal
from PyQt5.QtGui import QFont

REGION_THRESHOLDS = {
    "at_rich_upper": 47.5,
    "gc_rich_lower": 57.5,
    "neutral_lower": 42.5,
    "neutral_upper": 62.5,
}

PRIMER_DEFAULTS = {
    "AT-rich": {
        "outer": {
            "gc_min": 30,
            "gc_max": 65,
            "tm_min": 55.0,
            "tm_max": 58.0,
            "len_min": 18,
            "len_max": 25,
        },
        "inner": {
            "gc_min": 30,
            "gc_max": 65,
            "tm_min": 60.0,
            "tm_max": 63.0,
            "len_min": 20,
            "len_max": 25,
        },
        "loop": {
            "gc_min": 30,
            "gc_max": 65,
            "tm_min": 60.0,
            "tm_max": 63.0,
            "len_min": 20,
            "len_max": 25,
        },
    },
    "GC-rich": {
        "outer": {
            "gc_min": 40,
            "gc_max": 70,
            "tm_min": 59.0,
            "tm_max": 63.0,
            "len_min": 15,
            "len_max": 20,
        },
        "inner": {
            "gc_min": 40,
            "gc_max": 70,
            "tm_min": 64.0,
            "tm_max": 68.0,
            "len_min": 15,
            "len_max": 22,
        },
        "loop": {
            "gc_min": 40,
            "gc_max": 70,
            "tm_min": 64.0,
            "tm_max": 68.0,
            "len_min": 15,
            "len_max": 22,
        },
    },
    "Neutral": {
        "outer": {
            "gc_min": 40,
            "gc_max": 65,
            "tm_min": 59.0,
            "tm_max": 61.0,
            "len_min": 18,
            "len_max": 20,
        },
        "inner": {
            "gc_min": 40,
            "gc_max": 65,
            "tm_min": 64.0,
            "tm_max": 66.0,
            "len_min": 20,
            "len_max": 22,
        },
        "loop": {
            "gc_min": 40,
            "gc_max": 65,
            "tm_min": 64.0,
            "tm_max": 66.0,
            "len_min": 20,
            "len_max": 22,
        },
    },
}

REGION_HEADER_COLOR = {
    "AT-rich": "#dcdcdc",
    "GC-rich": "#dcdcdc",
    "Neutral": "#dcdcdc",
}
REGION_BODY_COLOR = {
    "AT-rich": "#e9e9e9",
    "GC-rich": "#e9e9e9",
    "Neutral": "#e9e9e9",
}


class RegionSection(QWidget):
    """
    Collapsible panel for one GC region (neutral, GC-rich, AT-rich)
    """

    changed = pyqtSignal()

    def __init__(self, region_name, defaults, always_enabled=False, parent=None):
        super().__init__(parent)
        self._expanded = False
        self.spinboxes = {"outer": {}, "inner": {}, "loop": {}}
        self._loop_widgets = []

        outer_lay = QVBoxLayout(self)
        outer_lay.setContentsMargins(0, 0, 0, 6)
        outer_lay.setSpacing(0)

        hc = REGION_HEADER_COLOR[region_name]
        bc = REGION_BODY_COLOR[region_name]

        header = QWidget()
        header.setStyleSheet(
            f"background:{hc}; border:1px solid #bbbbbb; border-radius:4px;"
        )
        h_lay = QHBoxLayout(header)
        h_lay.setContentsMargins(8, 6, 8, 6)
        h_lay.setSpacing(8)

        self._toggle_btn = QPushButton("▶")
        self._toggle_btn.setFixedSize(22, 22)
        self._toggle_btn.setFlat(True)
        self._toggle_btn.setStyleSheet("font-size:11px; color:#000000;")
        self._toggle_btn.clicked.connect(self._toggle)
        h_lay.addWidget(self._toggle_btn)

        name_lbl = QLabel(region_name)
        name_lbl.setStyleSheet("font-weight:bold; font-size:13px; color:#000000;")
        h_lay.addWidget(name_lbl)

        h_lay.addStretch()

        if always_enabled:
            tag = QLabel("Always on")
            tag.setStyleSheet("color:#333333; font-size:10px;")
            h_lay.addWidget(tag)
            self._enable_cb = None
        else:
            self._enable_cb = QCheckBox("Enabled")
            self._enable_cb.setStyleSheet("color:#000000;")
            self._enable_cb.setChecked(True)
            self._enable_cb.toggled.connect(self._on_enable_toggled)
            self._enable_cb.toggled.connect(self.changed.emit)
            h_lay.addWidget(self._enable_cb)

        outer_lay.addWidget(header)

        self._body = QWidget()
        self._body.setStyleSheet(
            f"background:{bc}; border:1px solid #bbbbbb; "
            "border-top:none; border-radius:0 0 4px 4px;"
        )
        grid = QGridLayout(self._body)
        grid.setContentsMargins(16, 10, 16, 10)
        grid.setHorizontalSpacing(14)
        grid.setVerticalSpacing(8)

        for col, hdr_text in enumerate(
            ["GC content (%)", "Tm (°C)", "Length (bp)"], start=1
        ):
            lbl = QLabel(hdr_text)
            lbl.setStyleSheet(
                "font-weight:bold; color:#333333; background-color:#ececec; "
                "border-radius:4px; padding:4px;"
            )
            lbl.setAlignment(Qt.AlignCenter)
            grid.addWidget(lbl, 0, col)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
        grid.addWidget(sep, 1, 0, 1, 4)

        row_defs = [
            ("Outer", "outer"),
            ("Inner", "inner"),
            ("Loop", "loop"),
        ]
        params = [
            ("gc_min", "gc_max", False),
            ("tm_min", "tm_max", True),
            ("len_min", "len_max", False),
        ]

        for grid_row, (label, ptype) in enumerate(row_defs, start=2):
            row_lbl = QLabel(label)
            row_lbl.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
            row_lbl.setStyleSheet(
                "font-size:12px; font-weight:bold; color:#333333; "
                "background-color:#ececec; border-radius:4px; padding:4px;"
            )
            grid.addWidget(row_lbl, grid_row, 0)

            for col, (k_min, k_max, is_float) in enumerate(params, start=1):
                d = defaults[ptype]
                sb_min = self._make_sb(d[k_min], is_float)
                sb_max = self._make_sb(d[k_max], is_float)
                cell = self._range_cell(sb_min, sb_max)
                grid.addWidget(cell, grid_row, col)

                self.spinboxes[ptype][k_min] = sb_min
                self.spinboxes[ptype][k_max] = sb_max
                sb_min.valueChanged.connect(self.changed.emit)
                sb_max.valueChanged.connect(self.changed.emit)

                if ptype == "loop":
                    self._loop_widgets.extend([sb_min, sb_max, row_lbl])

        self._body.setVisible(False)
        outer_lay.addWidget(self._body)

    @staticmethod
    def _make_sb(value, is_float):
        if is_float:
            sb = QDoubleSpinBox()
            sb.setRange(0.0, 200.0)
            sb.setSingleStep(0.5)
            sb.setDecimals(1)
            sb.setValue(float(value))
        else:
            sb = QSpinBox()
            sb.setRange(0, 200)
            sb.setValue(int(value))
        sb.setButtonSymbols(QAbstractSpinBox.NoButtons)
        sb.setStyleSheet("background-color:#ffffff;")
        sb.setFixedWidth(64)
        return sb

    @staticmethod
    def _range_cell(sb_min, sb_max):
        cell = QWidget()
        cell.setStyleSheet("background-color:#ececec; border-radius:4px;")
        lay = QHBoxLayout(cell)
        lay.setContentsMargins(6, 4, 6, 4)
        lay.setSpacing(4)
        lay.addStretch()
        lay.addWidget(sb_min)
        dash = QLabel("–")
        dash.setStyleSheet("background-color:transparent; border:none;")
        dash.setAlignment(Qt.AlignCenter)
        lay.addWidget(dash)
        lay.addWidget(sb_max)
        lay.addStretch()
        return cell

    def _toggle(self):
        self._expanded = not self._expanded
        self._body.setVisible(self._expanded)
        self._toggle_btn.setText("▼" if self._expanded else "▶")

    def _on_enable_toggled(self, checked):
        for ptype in ("outer", "inner", "loop"):
            for sb in self.spinboxes[ptype].values():
                sb.setEnabled(checked)

    def set_loop_enabled(self, enabled):
        """Called by the global loop checkbox to grey out / restore loop row."""
        region_enabled = self._enable_cb is None or self._enable_cb.isChecked()
        for w in self._loop_widgets:
            w.setEnabled(enabled and region_enabled)

    def get_values(self):
        return {
            ptype: {k: sb.value() for k, sb in sbs.items()}
            for ptype, sbs in self.spinboxes.items()
        }


class LAMPInterface(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("LAMP Primer Designer")
        self.setMinimumSize(820, 700)
        self.selected_file = None
        self.region_sections = {}
        self.threshold_boxes = {}
        self._region_buttons = []
        self._region_status_labels = []
        self._set_buttons = []
        self._set_status_labels = []
        self._proc = None
        self._proc_output_buffer = ""
        self._tmp_design = None
        self._candidates_fresh = False
        self._ever_ran = False
        self._resolved_input_path = None
        self._owns_input_path = False
        self._sets_proc = None
        self._sets_proc_output_buffer = ""
        self._results_out_path = None
        self._sets_fresh = False
        self._sets_ever_ran = False

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)

        self.tabs = QTabWidget()
        root.addWidget(self.tabs)
        self.tabs.addTab(self._make_sequence_tab(), "Sequence Input")
        self.tabs.addTab(self._make_candidate_params_tab(), "Candidate Region Params")
        self.tabs.addTab(self._make_set_params_tab(), "Set Params")
        self.tabs.addTab(self._make_results_tab(), "Results")
        self._refresh_set_buttons()

    # Sequence Input

    def _make_sequence_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setSpacing(10)

        file_group = QGroupBox("Input File  (.fa / .fna / .txt)")
        fg_lay = QHBoxLayout(file_group)

        self.file_label = QLabel("No file selected")
        self.file_label.setStyleSheet("color:grey; font-style:italic;")
        fg_lay.addWidget(self.file_label, stretch=1)

        browse_btn = QPushButton("Browse…")
        browse_btn.setFixedWidth(90)
        browse_btn.clicked.connect(self._browse_file)
        fg_lay.addWidget(browse_btn)

        clear_btn = QPushButton("Clear")
        clear_btn.setFixedWidth(70)
        clear_btn.clicked.connect(self._clear_file)
        fg_lay.addWidget(clear_btn)

        lay.addWidget(file_group)

        seq_group = QGroupBox("Or paste sequence directly (FASTA or raw nucleotides)")
        sg_lay = QVBoxLayout(seq_group)
        self.seq_text = QTextEdit()
        self.seq_text.setPlaceholderText(
            ">sequence_name\nACGTACGTACGT…\n\n" "Pasting here overrides selected file"
        )
        self.seq_text.setFont(QFont("Courier New", 10))
        self.seq_text.textChanged.connect(self._mark_dirty)
        sg_lay.addWidget(self.seq_text)
        lay.addWidget(seq_group, stretch=1)

        lay.addLayout(self._make_action_buttons_row())

        return w

    @staticmethod
    def _make_status_label():
        lbl = QLabel("")
        lbl.setAlignment(Qt.AlignCenter)
        lbl.setStyleSheet("font-size:11px; color:#555555;")
        lbl.setWordWrap(True)
        return lbl

    def _make_action_buttons_row(self):
        btn_lay = QHBoxLayout()
        btn_lay.addStretch()

        region_col = QVBoxLayout()
        region_btn = QPushButton("Create Candidate Regions")
        region_btn.setFixedWidth(200)
        region_btn.setFixedHeight(38)
        region_btn.setStyleSheet("font-weight:bold; font-size:13px;")
        region_btn.clicked.connect(self._run_candidates)
        self._region_buttons.append(region_btn)
        region_col.addWidget(region_btn)
        region_status = self._make_status_label()
        self._region_status_labels.append(region_status)
        region_col.addWidget(region_status)
        btn_lay.addLayout(region_col)

        btn_lay.addSpacing(20)

        sets_col = QVBoxLayout()
        sets_btn = QPushButton("Create Sets")
        sets_btn.setFixedWidth(160)
        sets_btn.setFixedHeight(38)
        sets_btn.setStyleSheet("font-weight:bold; font-size:13px;")
        sets_btn.clicked.connect(self._run_sets)
        self._set_buttons.append(sets_btn)
        sets_col.addWidget(sets_btn)
        sets_status = self._make_status_label()
        self._set_status_labels.append(sets_status)
        sets_col.addWidget(sets_status)
        btn_lay.addLayout(sets_col)

        btn_lay.addStretch()
        return btn_lay

    def _browse_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select sequence file",
            "",
            "Sequence files (*.fa *.fna *.fasta *.txt);;All files (*)",
        )
        if path:
            self.selected_file = path
            self.file_label.setText(os.path.basename(path))
            self.file_label.setStyleSheet("color:black; font-style:normal;")
            self._mark_dirty()

    def _clear_file(self):
        self.selected_file = None
        self.file_label.setText("No file selected")
        self.file_label.setStyleSheet("color:grey; font-style:italic;")
        self._mark_dirty()

    def _run_candidates(self):
        if self._proc is not None:
            return

        pasted = self.seq_text.toPlainText().strip()
        has_file = bool(self.selected_file)
        if not has_file and not pasted:
            self.file_label.setText(
                "Please select a sequence file or paste a sequence."
            )
            self.file_label.setStyleSheet("color:red; font-style:italic;")
            self.tabs.setCurrentIndex(0)
            return

        self._discard_resolved_input()
        if pasted:
            text = pasted if pasted.startswith(">") else ">pasted_sequence\n" + pasted
            tmp_fa = tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False
            )
            tmp_fa.write(text)
            tmp_fa.close()
            self._resolved_input_path = tmp_fa.name
            self._owns_input_path = True
        else:
            self._resolved_input_path = self.selected_file
            self._owns_input_path = False

        tmp_design = tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False)
        json.dump(self._build_design(), tmp_design, indent=2)
        tmp_design.close()
        self._tmp_design = tmp_design.name

        script_dir = os.path.dirname(os.path.abspath(__file__))
        script = os.path.join(script_dir, "single.py")
        args = [
            script,
            "-in",
            self._resolved_input_path,
            "-design",
            self._tmp_design,
            "-loop",
            "True" if self.loop_cb.isChecked() else "False",
        ]

        self._proc_output_buffer = ""
        self._ever_ran = True
        self._proc = QProcess(self)
        self._refresh_region_buttons()

        self._proc.setWorkingDirectory(script_dir)
        self._proc.setProcessChannelMode(QProcess.MergedChannels)
        self._proc.readyReadStandardOutput.connect(self._on_proc_output)
        self._proc.finished.connect(self._on_proc_finished)
        self._proc.errorOccurred.connect(self._on_proc_error)
        self._proc.start(sys.executable, args)

    def _build_design(self):
        regions = {
            region: sec.get_values() for region, sec in self.region_sections.items()
        }
        gc_auto = self.gc_auto_cb.isChecked()
        return {
            "region_thresholds": {
                k: sb.value() for k, sb in self.threshold_boxes.items()
            },
            "regions": regions,
            "global": {
                "min_len": self.len_min_sb.value(),
                "max_len": self.len_max_sb.value(),
                "tm_min": self.tm_min_sb.value(),
                "tm_max": self.tm_max_sb.value(),
                "gc_min": None if gc_auto else self.gc_min_sb.value(),
                "gc_max": None if gc_auto else self.gc_max_sb.value(),
            },
        }

    def _on_proc_output(self):
        self._proc_output_buffer += bytes(self._proc.readAllStandardOutput()).decode(
            errors="replace"
        )

    def _set_region_status(self, text, color):
        for lbl in self._region_status_labels:
            lbl.setText(text)
            lbl.setStyleSheet(f"font-size:11px; color:{color};")

    def _on_proc_finished(self, exit_code, _exit_status):
        if exit_code == 0:
            self._candidates_fresh = True
            self._sets_fresh = False
            summary = self._format_candidate_summary(self._proc_output_buffer)
            self._set_region_status(summary or "No candidates found in output", "green")
        else:
            self._candidates_fresh = False
            tail = self._proc_output_buffer.strip().splitlines()[-1:] or [""]
            self._set_region_status(f"Failed", "darkred")

        self._cleanup_run()

    def _on_proc_error(self, error):
        self._candidates_fresh = False
        self._set_region_status(f"Process error: {error}", "darkred")
        self._cleanup_run()

    def _format_candidate_summary(self, text):
        m = re.search(
            r"Total candidates:\s*(\d+)\s*inner,\s*(\d+)\s*outer(?:,\s*(\d+)\s*loop)?",
            text,
        )
        if not m:
            return None
        inner, outer, loop = m.groups()
        parts = f"{inner} inner, {outer} outer"
        if loop is not None:
            parts += f", {loop} loop"

        done = re.search(r"Done in ([\d.]+)s", text)
        time_str = f"{done.group(1)}s" if done else "?"
        return f"{parts} candidates created in {time_str}"

    def _discard_resolved_input(self):
        if (
            self._owns_input_path
            and self._resolved_input_path
            and os.path.exists(self._resolved_input_path)
        ):
            os.remove(self._resolved_input_path)
        self._resolved_input_path = None
        self._owns_input_path = False

    def _cleanup_run(self):
        self._proc = None
        self._refresh_region_buttons()
        self._refresh_set_buttons()
        if self._tmp_design and os.path.exists(self._tmp_design):
            os.remove(self._tmp_design)
        self._tmp_design = None

    def _mark_dirty(self, *_args):
        was_fresh = self._candidates_fresh
        self._candidates_fresh = False
        self._sets_fresh = False
        if (was_fresh or self._ever_ran) and self._region_status_labels:
            self._set_region_status("Parameters changed, rerun", "sienna")
        self._refresh_region_buttons()
        self._refresh_set_buttons()

    def _mark_sets_dirty(self, *_args):
        was_fresh = self._sets_fresh
        self._sets_fresh = False
        if (was_fresh or self._sets_ever_ran) and self._set_status_labels:
            self._set_sets_status("Parameters changed, rerun", "sienna")
        self._refresh_set_buttons()

    def _refresh_region_buttons(self):
        running = self._proc is not None or self._sets_proc is not None
        fresh = self._candidates_fresh
        for btn in self._region_buttons:
            if running:
                btn.setEnabled(False)
                btn.setText("Running…")
            elif fresh:
                btn.setEnabled(False)
                btn.setText("Candidate Regions Ready")
            else:
                btn.setEnabled(True)
                btn.setText("Create Candidate Regions")

    def _run_sets(self):
        if self._sets_proc is not None or self._proc is not None:
            return
        if not self._candidates_fresh or not self._resolved_input_path:
            return

        script_dir = os.path.dirname(os.path.abspath(__file__))
        script = os.path.join(script_dir, "lamp.py")
        args = [
            script,
            "-in",
            self._resolved_input_path,
            "-num",
            str(self.num_sb.value()),
            "-gblock",
            "True" if self.gblock_cb.isChecked() else "False",
            "-loop",
            "True" if self.loop_cb.isChecked() else "False",
            "-multiplex",
            str(self.multiplex_sb.value()),
            "-max_repairs",
            str(self.max_repairs_sb.value()),
        ]

        self._sets_proc_output_buffer = ""
        self._sets_ever_ran = True
        self._sets_proc = QProcess(self)
        self._refresh_region_buttons()
        self._refresh_set_buttons()

        self._sets_proc.setWorkingDirectory(script_dir)
        self._sets_proc.setProcessChannelMode(QProcess.MergedChannels)
        self._sets_proc.readyReadStandardOutput.connect(self._on_sets_proc_output)
        self._sets_proc.finished.connect(self._on_sets_proc_finished)
        self._sets_proc.errorOccurred.connect(self._on_sets_proc_error)
        self._sets_proc.start(sys.executable, args)

    def _on_sets_proc_output(self):
        self._sets_proc_output_buffer += bytes(
            self._sets_proc.readAllStandardOutput()
        ).decode(errors="replace")

    def _set_sets_status(self, text, color):
        for lbl in self._set_status_labels:
            lbl.setText(text)
            lbl.setStyleSheet(f"font-size:11px; color:{color};")

    def _on_sets_proc_finished(self, exit_code, _exit_status):
        if exit_code == 0:
            self._sets_fresh = True
            summary = self._format_sets_summary(self._sets_proc_output_buffer)
            self._set_sets_status(summary or "No sets found", "green")
            script_dir = os.path.dirname(os.path.abspath(__file__))
            results_path = os.path.join(script_dir, "Results", "results.txt")
            if os.path.exists(results_path):
                with open(results_path) as fh:
                    self.results_text.setPlainText(fh.read())
            self.save_results_btn.setEnabled(
                bool(self.results_text.toPlainText().strip())
            )
            if self._results_out_path:
                self._save_results(self._results_out_path)
        else:
            self._sets_fresh = False
            tail = self._sets_proc_output_buffer.strip().splitlines()[-1:] or [""]
            self._set_sets_status(
                f"Failed (exit code {exit_code}): {tail[0]}"[:140], "darkred"
            )

        self._cleanup_sets_run()

    def _on_sets_proc_error(self, error):
        self._sets_fresh = False
        self._set_sets_status(f"Process error: {error}", "darkred")
        self._cleanup_sets_run()

    def _format_sets_summary(self, text):
        m = re.search(r"Assembly complete:\s*(\d+)\s*sets found in ([\d.]+)s", text)
        if not m:
            return None
        n, time_str = m.groups()
        return f"{n} sets created in {time_str}s"

    def _cleanup_sets_run(self):
        self._sets_proc = None
        self._refresh_region_buttons()
        self._refresh_set_buttons()

    def _refresh_set_buttons(self):
        running = self._sets_proc is not None
        blocked = self._proc is not None or not self._candidates_fresh
        for btn in self._set_buttons:
            if running:
                btn.setEnabled(False)
                btn.setText("Running…")
            elif self._sets_fresh:
                btn.setEnabled(False)
                btn.setText("Sets Ready")
            elif blocked:
                btn.setEnabled(False)
                btn.setText("Create Sets")
            else:
                btn.setEnabled(True)
                btn.setText("Create Sets")

    # Candidate Regions Parameters

    def _make_candidate_params_tab(self):
        outer_w = QWidget()
        outer_lay = QVBoxLayout(outer_w)
        outer_lay.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer_lay.addWidget(scroll)

        container = QWidget()
        scroll.setWidget(container)
        c_lay = QVBoxLayout(container)
        c_lay.setSpacing(12)
        c_lay.setContentsMargins(10, 10, 10, 10)

        c_lay.addWidget(self._make_threshold_group())
        c_lay.addWidget(self._make_global_group())

        # Collapsible region sections

        for region in ("Neutral", "AT-rich", "GC-rich"):
            sec = RegionSection(
                region,
                PRIMER_DEFAULTS[region],
                always_enabled=(region == "Neutral"),
            )
            self.region_sections[region] = sec
            sec.changed.connect(self._mark_dirty)
            c_lay.addWidget(sec)

        self.loop_cb.toggled.connect(self._on_loop_toggled)

        c_lay.addStretch()

        btn_lay = self._make_action_buttons_row()
        btn_lay.setContentsMargins(10, 6, 10, 6)
        outer_lay.addLayout(btn_lay)

        return outer_w

    def _make_threshold_group(self):
        grp = QGroupBox("GC Region Thresholds")
        lay = QVBoxLayout(grp)
        lay.setSpacing(10)

        thresh_row = QHBoxLayout()
        thresh_row.setSpacing(10)

        thresh_row.addWidget(QLabel("AT-rich if regional GC% <"))
        self.at_sb = self._make_float_sb(
            REGION_THRESHOLDS["at_rich_upper"], suffix=" %"
        )
        thresh_row.addWidget(self.at_sb)

        thresh_row.addSpacing(20)

        thresh_row.addWidget(QLabel("GC-rich if regional GC% >"))
        self.gc_sb = self._make_float_sb(
            REGION_THRESHOLDS["gc_rich_lower"], suffix=" %"
        )
        thresh_row.addWidget(self.gc_sb)

        thresh_row.addStretch()
        lay.addLayout(thresh_row)

        neutral_row = QHBoxLayout()
        neutral_row.setSpacing(10)

        neutral_row.addWidget(QLabel("Neutral range:"))
        self.neutral_lo_sb = self._make_float_sb(
            REGION_THRESHOLDS["neutral_lower"], suffix=" %"
        )
        neutral_row.addWidget(self.neutral_lo_sb)
        neutral_row.addWidget(QLabel("–"))
        self.neutral_hi_sb = self._make_float_sb(
            REGION_THRESHOLDS["neutral_upper"], suffix=" %"
        )
        neutral_row.addWidget(self.neutral_hi_sb)

        neutral_row.addSpacing(10)
        neutral_row.addStretch()
        lay.addLayout(neutral_row)

        self.threshold_boxes = {
            "at_rich_upper": self.at_sb,
            "gc_rich_lower": self.gc_sb,
            "neutral_lower": self.neutral_lo_sb,
            "neutral_upper": self.neutral_hi_sb,
        }
        for sb in self.threshold_boxes.values():
            sb.valueChanged.connect(self._mark_dirty)
        return grp

    def _make_global_group(self):
        grp = QGroupBox("Global Pre-filter")
        lay = QVBoxLayout(grp)
        lay.setSpacing(10)

        row1 = QHBoxLayout()
        row1.setSpacing(10)
        row1.addWidget(QLabel("Primer length:"))
        self.len_min_sb = RegionSection._make_sb(15, False)
        row1.addWidget(self.len_min_sb)
        row1.addWidget(QLabel("–"))
        self.len_max_sb = RegionSection._make_sb(25, False)
        row1.addWidget(self.len_max_sb)
        row1.addWidget(QLabel("bp"))

        row1.addSpacing(20)
        row1.addWidget(QLabel("Tm:"))
        self.tm_min_sb = self._make_float_sb(55.0)
        row1.addWidget(self.tm_min_sb)
        row1.addWidget(QLabel("–"))
        self.tm_max_sb = self._make_float_sb(68.0)
        row1.addWidget(self.tm_max_sb)
        row1.addWidget(QLabel("°C"))
        row1.addStretch()
        lay.addLayout(row1)

        row2 = QHBoxLayout()
        row2.setSpacing(10)
        row2.addWidget(QLabel("GC content:"))
        self.gc_min_sb = RegionSection._make_sb(30, False)
        row2.addWidget(self.gc_min_sb)
        row2.addWidget(QLabel("–"))
        self.gc_max_sb = RegionSection._make_sb(70, False)
        row2.addWidget(self.gc_max_sb)
        row2.addWidget(QLabel("%"))

        self.gc_auto_cb = QCheckBox("Scale with sequence length")
        self.gc_auto_cb.setChecked(True)
        self.gc_auto_cb.toggled.connect(self._on_gc_auto_toggled)
        row2.addSpacing(15)
        row2.addWidget(self.gc_auto_cb)
        row2.addStretch()
        lay.addLayout(row2)

        loop_row = QHBoxLayout()
        self.loop_cb = QCheckBox("Design loop primers (LF/LB)")
        self.loop_cb.setChecked(True)
        loop_row.addWidget(self.loop_cb)
        loop_row.addStretch()
        lay.addLayout(loop_row)

        self.gc_min_sb.setEnabled(False)
        self.gc_max_sb.setEnabled(False)

        for sb in (
            self.len_min_sb,
            self.len_max_sb,
            self.tm_min_sb,
            self.tm_max_sb,
            self.gc_min_sb,
            self.gc_max_sb,
        ):
            sb.valueChanged.connect(self._mark_dirty)
        self.gc_auto_cb.toggled.connect(self._mark_dirty)
        self.loop_cb.toggled.connect(self._mark_dirty)

        return grp

    def _on_gc_auto_toggled(self, checked):
        self.gc_min_sb.setEnabled(not checked)
        self.gc_max_sb.setEnabled(not checked)

    def _on_loop_toggled(self, enabled):
        for sec in self.region_sections.values():
            sec.set_loop_enabled(enabled)

    def _make_set_params_tab(self):
        outer_w = QWidget()
        outer_lay = QVBoxLayout(outer_w)
        outer_lay.setContentsMargins(0, 0, 0, 0)

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        outer_lay.addWidget(scroll)

        container = QWidget()
        scroll.setWidget(container)
        c_lay = QVBoxLayout(container)
        c_lay.setSpacing(12)
        c_lay.setContentsMargins(10, 10, 10, 10)

        c_lay.addWidget(self._make_set_count_group())
        c_lay.addWidget(self._make_repairs_group())
        c_lay.addWidget(self._make_gblock_group())
        c_lay.addStretch()

        btn_lay = self._make_action_buttons_row()
        btn_lay.setContentsMargins(10, 6, 10, 6)
        outer_lay.addLayout(btn_lay)

        return outer_w

    def _make_set_count_group(self):
        grp = QGroupBox("Set Count")
        lay = QHBoxLayout(grp)
        lay.addWidget(QLabel("Number of LAMP sets to find:"))
        self.num_sb = QSpinBox()
        self.num_sb.setRange(1, 1000)
        self.num_sb.setValue(10)
        self.num_sb.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.num_sb.setFixedWidth(70)
        self.num_sb.valueChanged.connect(self._mark_sets_dirty)
        lay.addWidget(self.num_sb)
        lay.addStretch()
        return grp

    def _make_repairs_group(self):
        grp = QGroupBox("Repairs")
        lay = QVBoxLayout(grp)
        lay.setSpacing(10)

        rep_row = QHBoxLayout()
        rep_row.setSpacing(10)
        rep_row.addWidget(QLabel("Allowed degenerate bases (0 for no repair):"))
        self.max_repairs_sb = QSpinBox()
        self.max_repairs_sb.setRange(0, 20)
        self.max_repairs_sb.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.max_repairs_sb.setFixedWidth(70)
        self.max_repairs_sb.valueChanged.connect(self._mark_sets_dirty)
        rep_row.addWidget(self.max_repairs_sb)
        rep_row.addStretch()
        lay.addLayout(rep_row)

        ms_row = QHBoxLayout()
        ms_row.setSpacing(10)
        ms_row.addWidget(QLabel("Multiplex bundle size (0 for single sets):"))
        self.multiplex_sb = QSpinBox()
        self.multiplex_sb.setRange(0, 20)
        self.multiplex_sb.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.multiplex_sb.setFixedWidth(70)
        self.multiplex_sb.valueChanged.connect(self._mark_sets_dirty)
        ms_row.addWidget(self.multiplex_sb)
        ms_row.addStretch()
        lay.addLayout(ms_row)

        return grp

    def _make_gblock_group(self):
        grp = QGroupBox("gBlock")
        lay = QHBoxLayout(grp)
        self.gblock_cb = QCheckBox("Design/require a synthesizable gBlock for each set")
        self.gblock_cb.setChecked(True)
        self.gblock_cb.toggled.connect(self._mark_sets_dirty)
        lay.addWidget(self.gblock_cb)
        lay.addStretch()
        return grp

    def _make_results_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(10, 10, 10, 10)
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setFont(QFont("Courier New", 10))
        self.results_text.setPlaceholderText(
            "Sets will appear here after Create Sets is run"
        )
        lay.addWidget(self.results_text)

        save_row = QHBoxLayout()
        self.save_results_btn = QPushButton("Save Results As…")
        self.save_results_btn.setEnabled(False)
        self.save_results_btn.clicked.connect(self._choose_results_output)
        self.results_out_label = QLabel("Not saved")
        self.results_out_label.setStyleSheet("font-size:11px; color:#666666;")
        save_row.addWidget(self.save_results_btn)
        save_row.addWidget(self.results_out_label, 1)
        lay.addLayout(save_row)
        return w

    def _save_results(self, path):
        text = self.results_text.toPlainText()
        if not text.strip():
            return False
        parent = os.path.dirname(path) or "."
        if not os.access(parent, os.W_OK):
            self.results_out_label.setText(f"Save failed", "darkred")
            self.results_out_label.setStyleSheet("font-size:11px; color:darkred;")
            return False
        with open(path, "w") as fh:
            fh.write(text)
        self._results_out_path = path
        self.results_out_label.setText(f"Saved to {path}")
        self.results_out_label.setStyleSheet("font-size:11px; color:green;")
        return True

    def _choose_results_output(self):
        if not self.results_text.toPlainText().strip():
            return
        start = self._results_out_path or os.path.join(
            os.path.expanduser("~"), "lamp_results.txt"
        )
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Results", start, "Text files (*.txt);;All files (*)"
        )
        if path:
            self._save_results(path)

    @staticmethod
    def _make_float_sb(value, suffix="", width=80):
        sb = QDoubleSpinBox()
        sb.setRange(0.0, 100.0)
        sb.setSingleStep(0.5)
        sb.setDecimals(1)
        sb.setValue(float(value))
        if suffix:
            sb.setSuffix(suffix)
        sb.setButtonSymbols(QAbstractSpinBox.NoButtons)
        sb.setFixedWidth(width)
        return sb


def main():
    app = QApplication(sys.argv)
    win = LAMPInterface()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
