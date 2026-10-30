#!/usr/bin/env python3
"""Small integration tests for arbitrary image counts and unchanged core methods.

Run from the repository root: python tests/test_n_images.py
The fixtures are distinct synthetic 2D fields, not independent scientific
validation of any material. All generated files live in a temporary directory.
"""
from pathlib import Path
import csv
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from PIL import Image

REPO = Path(__file__).resolve().parents[1]


class NImageIntegration(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory(prefix="rea-n-images-")
        cls.work = Path(cls.tmp.name)
        cls.inputs = cls.work / "input fields"
        cls.inputs.mkdir()
        cls.images, cls.masks, cls.kappa, cls.young, cls.poisson = [], [], [], [], []
        rng = np.random.default_rng(4056)
        for i in range(20):
            # Distinct values and image names deliberately unrelated to 1_cut, etc.
            stem = f"sample_{chr(65 + i)}"
            image = cls.inputs / f"{stem}.tif"
            gray = rng.integers(65 + i, 155 + i, size=(96, 96), dtype=np.uint8)
            Image.fromarray(gray).save(image)
            mask = cls.inputs / f"{stem}_stationary_mask.tif"
            Image.fromarray(np.full(gray.shape, 255, dtype=np.uint8)).save(mask)
            cls.images.append(image)
            cls.masks.append(mask)
            # Homogeneous cells have an exact thermal solution and a known
            # plane-strain axial compliance, so the FEM test checks values.
            for group, prefix, value in [(cls.kappa, "conductivity", 1.0 + .05*i),
                                         (cls.young, "elastic", 10.0 + i),
                                         (cls.poisson, "poisson", .20 + .002*i)]:
                p = cls.inputs / f"{prefix}_{stem}.csv"
                np.savetxt(p, np.full((16, 16), value), delimiter=",")
                group.append(p)
        cls.empty_mask = cls.inputs / "empty_mask.tif"
        Image.fromarray(np.zeros((96, 96), dtype=np.uint8)).save(cls.empty_mask)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def run_program(self, script, *args, succeeds=True):
        env = os.environ.copy()
        env.pop("REA_DATA_ROOT", None)
        env["MPLBACKEND"] = "Agg"
        result = subprocess.run(
            [sys.executable, str(REPO / script), *map(str, args)],
            cwd=self.work, env=env, capture_output=True, text=True, timeout=180,
        )
        if succeeds:
            self.assertEqual(result.returncode, 0, result.stdout[-4000:] + result.stderr[-4000:])
        else:
            self.assertNotEqual(result.returncode, 0, "Invalid input was silently accepted")
        return result

    @staticmethod
    def read_csv(path):
        with Path(path).open(newline="") as handle:
            return list(csv.DictReader(handle))

    def validation_args(self, n, output, masks=None):
        return ["--images", *self.images[:n], "--masks", *(masks or self.masks[:n]),
                "--window-sizes", "8,16,24,32", "--reference-L", "32",
                "--extra-reference-sizes", "48,64", "--n-window-samples", "9",
                "--n-common-k-grid", "32", "--tau-Z", "2", "--tau-C", "10",
                "--screening-sizes", "16,24", "--screening-stride", "8",
                "--screening-tau-mu", "2", "--screening-tau-sigma", "2",
                "--min-component-windows", "1",
                "--output-dir", output, "--checkpoint-dir", output / "checkpoints"]

    def test_core_sources_unchanged(self):
        expected = {
            "detect_stationary_domains_general.py": "9682ff8bec2d55e2f2481fb361326ed47ee8c943",
            "rea_stationary_mask_cli_v24_general.py": "862983356fc38aa8df56a1a3dea1fc962c150537",
            "property_full_homogenized_inside_property_masks.py": "ca4ae44afa2e691255608a4d62e38dbd3ef92997",
        }
        for name, digest in expected.items():
            data = (REPO/name).read_bytes()
            self.assertEqual(hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest(), digest)

    def test_core_nine_and_twenty(self):
        for n in (9, 20):
            with self.subTest(n=n):
                out = self.work / f"core_{n}"
                detector = out / "masks"
                self.run_program(
                    "detect_stationary_domains_general.py", *self.images[:n],
                    "--input-type", "image", "--output-dir", detector,
                    "--stationarity-window", "16", "--stride", "8",
                    "--tau-mu", "2", "--tau-sigma", "2",
                    "--min-component-windows", "1", "--no-morphology", "--dpi", "35",
                )
                masks = sorted(detector.glob("*_stationary_mask.tif"))
                self.assertEqual(len(masks), n)
                for m in masks:
                    self.assertTrue(np.asarray(Image.open(m)).all())
                structural = out / "structural"
                self.run_program(
                    "rea_stationary_mask_cli_v24_general.py", *self.images[:n],
                    "--max-side", "0", "--window-sizes", "8,16,24,32",
                    "--reference-L", "32", "--selected-spectra", "8,32",
                    "--n-common-k-grid", "32", "--n-window-samples", "9",
                    "--tau-Z", "2", "--tau-C", "10", "--output-dir", structural,
                    "--output-prefix", "synthetic", "--stationary-mask-dir", detector,
                    "--dpi", "35", "--no-pdf",
                )
                rows = self.read_csv(structural / "synthetic_combined_window_statistics.csv")
                self.assertEqual(len(rows), 4)
                self.assertTrue(all(int(float(r["n_images"])) == n for r in rows))
                self.assertEqual(sum(k.startswith("sample_mean_Z__") for k in rows[0]), n)
                prop_masks = out / "property_masks"
                self.run_program(
                    "detect_stationary_domains_general.py", "--kappa-maps", *self.kappa[:n],
                    "--E-maps", *self.young[:n], "--nu-maps", *self.poisson[:n],
                    "--stationarity-source", "combined", "--output-dir", prop_masks,
                    "--stationarity-window", "8", "--stride", "4", "--tau-mu", "2",
                    "--tau-sigma", "2", "--min-component-windows", "1", "--no-morphology", "--dpi", "35",
                )
                props = out / "properties"
                self.run_program(
                    "property_full_homogenized_inside_property_masks.py", "--kappa-maps", *self.kappa[:n],
                    "--E-maps", *self.young[:n], "--nu-maps", *self.poisson[:n],
                    "--mask-dir", prop_masks, "--window-sizes", "4,8", "--L-rea", "4",
                    "--n-window-samples", "1", "--output-dir", props, "--output-prefix", "synthetic",
                    "--elastic-mode", "plane_strain", "--dpi", "35", "--no-pdf",
                )
                rows = self.read_csv(props / "synthetic_window_values.csv")
                self.assertEqual(len(rows), 2*n)
                self.assertEqual({int(r["image"]) for r in rows}, set(range(1, n+1)))
                for r in rows:
                    i = int(r["image"]) - 1
                    expected_k = 1.0 + .05*i
                    expected_E = (10.0+i)/(1-(.20+.002*i)**2)
                    np.testing.assert_allclose([float(r["K11"]),float(r["K22"])], expected_k, rtol=1e-10)
                    np.testing.assert_allclose([float(r["Ex"]),float(r["Ey"])], expected_E, rtol=1e-10)

    def test_validation_nine_and_twenty(self):
        for n in (9, 20):
            with self.subTest(n=n):
                out = self.work / f"validation_{n}"
                args = self.validation_args(n, out)
                self.run_program("validation/rea_validation_wrapper.py", "--all-images", *args)
                self.run_program("validation/rea_validation_wrapper.py", "--aggregate", "--bootstrap", "30",
                                 "--checkpoint-dir", out / "checkpoints", "--output-dir", out)
                result = json.loads((out / "validation_results.json").read_text())
                self.assertEqual(len(result["mask_fractions"]), n)
                self.assertEqual(len(result["leave_one_image_out"]), n)
                self.assertEqual(result["baseline"]["reference"], 32)
                self.assertEqual(result["baseline"]["sizes"], [8,16,24,32])
                with np.load(out / "bootstrap_draws.npz", allow_pickle=False) as draws:
                    self.assertEqual(draws["draws"].shape, (30, n))
                    self.assertEqual(set(draws["draws"].ravel()), set(range(n)))
                    self.assertEqual(draws["selected_pixels"].shape, (30,))
                self.run_program("validation/rebuild_derived.py", "--checkpoint-dir", out / "checkpoints", "--output-dir", out)
                with np.load(out / "figure2_plot_data.npz", allow_pickle=False) as plot:
                    np.testing.assert_array_equal(plot["L"], [8,16,24,32])
                    self.assertEqual(int(plot["reference_L"]), 32)
                    self.assertNotEqual(plot["accepted_L"].dtype, np.dtype("O"))
                self.run_program("validation/rstat_sensitivity_wrapper.py", "--all-images", *args)
                self.run_program("validation/rstat_sensitivity_wrapper.py", "--aggregate",
                                 "--checkpoint-dir", out / "checkpoints", "--output-dir", out)
                sensitivity = json.loads((out / "rstat_sensitivity_results.json").read_text())
                self.assertEqual(len(sensitivity["elapsed_image_seconds"]), n)
                self.assertEqual([r["R"] for r in sensitivity["experiments"]], [16,24])
                for experiment in sensitivity["experiments"]:
                    self.assertEqual(len(experiment["mask_fractions"]), n)
                # Reusing a set for different images must not mix old statistics.
                self.run_program("validation/rea_validation_wrapper.py", "--all-images",
                                 *self.validation_args(n-1, out), succeeds=False)
                changed = self.validation_args(n, out)
                changed[changed.index("--reference-L")+1] = "24"
                self.run_program("validation/rea_validation_wrapper.py", "--all-images", *changed, succeeds=False)
                # A missing member of an explicitly registered N-image dataset
                # must not silently turn its ensemble into N-1 images.
                checkpoints = sorted((out/"checkpoints").glob("image_*.json.gz"))
                self.assertEqual(len(checkpoints), n)
                if n == 20:
                    # Even the same image ID from another dataset is not an
                    # interchangeable checkpoint: its ensemble identity differs.
                    first = out/"checkpoints/image_1.json.gz"
                    original = first.read_bytes()
                    first.write_bytes((self.work/"validation_9/checkpoints/image_1.json.gz").read_bytes())
                    self.run_program("validation/rea_validation_wrapper.py", "--aggregate", "--bootstrap", "10",
                                     "--checkpoint-dir", out/"checkpoints", "--output-dir", out, succeeds=False)
                    first.write_bytes(original)
                checkpoints[-1].unlink()
                self.run_program("validation/rea_validation_wrapper.py", "--aggregate", "--bootstrap", "10",
                                 "--checkpoint-dir", out/"checkpoints", "--output-dir", out, succeeds=False)

    def test_no_selection_and_single_image(self):
        out = self.work / "none_selected"
        self.run_program("validation/rea_validation_wrapper.py", "--all-images",
                         *self.validation_args(9, out, [self.empty_mask]*9))
        self.run_program("validation/rea_validation_wrapper.py", "--aggregate", "--bootstrap", "20",
                         "--checkpoint-dir", out / "checkpoints", "--output-dir", out)
        result = json.loads((out / "validation_results.json").read_text())
        self.assertIsNone(result["baseline"]["selected"])
        self.run_program("validation/rebuild_derived.py", "--checkpoint-dir", out / "checkpoints", "--output-dir", out)
        with np.load(out / "figure2_plot_data.npz", allow_pickle=False) as plot:
            self.assertTrue(np.isnan(plot["accepted_L"]))
        self.run_program("figures/regenerate_figure2_layout.py", "--image-file", self.images[0],
                         "--data-file", out/"figure2_plot_data.npz", "--output-dir", out/"figure", "--dpi", "35")
        audit = json.loads((out/"figure/Figure_2_layout_audit.json").read_text())
        self.assertIsNone(audit["accepted_L"])
        one = self.work / "single"
        self.run_program("validation/rea_validation_wrapper.py", "--all-images", *self.validation_args(1, one))
        self.run_program("validation/rea_validation_wrapper.py", "--aggregate", "--bootstrap", "10",
                         "--checkpoint-dir", one / "checkpoints", "--output-dir", one)
        result = json.loads((one / "validation_results.json").read_text())
        self.assertEqual(result["leave_one_image_out"], [])
        self.assertFalse(result["bootstrap"]["independent_image_uncertainty_estimable"])
        self.assertIn("degenerate", result["bootstrap"]["note"])
        with np.load(one / "bootstrap_draws.npz", allow_pickle=False) as draws:
            self.assertEqual(draws["draws"].shape, (10, 1))

    def test_checkpoint_file_edits_rejected(self):
        for changed_input in ("image", "mask"):
            with self.subTest(changed_input=changed_input):
                out = self.work / f"changed_{changed_input}"
                out.mkdir()
                image, mask = out/"image.tif", out/"mask.tif"
                image.write_bytes(self.images[0].read_bytes())
                mask.write_bytes(self.masks[0].read_bytes())
                args = self.validation_args(1, out)
                args[args.index("--images")+1] = image
                args[args.index("--masks")+1] = mask
                self.run_program("validation/rea_validation_wrapper.py", "--all-images", *args)
                changed = image if changed_input == "image" else mask
                array = np.array(Image.open(changed))
                array[0, 0] = 0
                Image.fromarray(array).save(changed)
                self.run_program("validation/rea_validation_wrapper.py", "--all-images", *args, succeeds=False)

    def test_constant_image_checkpoint_roundtrip(self):
        image = self.inputs/"constant.tif"
        Image.fromarray(np.full((96, 96), 90, dtype=np.uint8)).save(image)
        for empty in (False, True):
            with self.subTest(empty_mask=empty):
                out = self.work / f"constant_{empty}"
                args = self.validation_args(1, out, [self.empty_mask] if empty else self.masks[:1])
                args[args.index("--images")+1] = image
                self.run_program("validation/rea_validation_wrapper.py", "--all-images", *args)
                # Constant fields create undefined boundary/variance ratios.
                # Those NaNs must survive the safe checkpoint format and be
                # consumable by aggregation and the plotting-data rebuild.
                self.run_program("validation/rea_validation_wrapper.py", "--aggregate", "--bootstrap", "10",
                                 "--checkpoint-dir", out/"checkpoints", "--output-dir", out)
                result = json.loads((out/"validation_results.json").read_text(),
                                    parse_constant=lambda value: self.fail("Non-standard JSON constant: " + value))
                self.assertFalse(result["bootstrap"]["independent_image_uncertainty_estimable"])
                self.run_program("validation/rebuild_derived.py", "--checkpoint-dir", out/"checkpoints", "--output-dir", out)
                with np.load(out/"figure2_plot_data.npz", allow_pickle=False) as plot:
                    np.testing.assert_array_equal(plot["L"], [8, 16, 24, 32])
                    self.assertNotEqual(plot["accepted_L"].dtype, np.dtype("O"))
                    if empty:
                        self.assertTrue(np.isnan(plot["accepted_L"]))

    def test_generic_structural_figures(self):
        out = self.work/"figures"
        image = self.inputs/"portrait.tif"
        Image.fromarray(np.arange(96*48,dtype=np.uint16).reshape(96,48).astype(np.uint8)).save(image)
        detector = out/"detector"
        self.run_program("detect_stationary_domains_general.py", image,
                         "--stationarity-window", "16", "--stride", "8", "--tau-mu", "2", "--tau-sigma", "2",
                         "--min-component-windows", "1", "--no-morphology", "--output-dir", detector, "--dpi", "35")
        self.run_program("figures/regenerate_figure1_layout.py", "--image-file", image,
                         "--mask-file", detector/"portrait_stationary_mask.tif",
                         "--windows-file", detector/"portrait_stationary_windows.csv",
                         "--output-dir", out/"fig1", "--dpi", "35")
        audit = json.loads((out/"fig1/Figure_1_layout_audit.json").read_text())
        self.assertEqual(audit["image_shape"], [96,48])
        self.assertTrue(audit["all_panel_dimensions_equal"])
        # A non-study support grid, arbitrary spectra, and an unselected support
        # check input-driven plot behavior without inventing a study result.
        k = np.linspace(0, np.pi, 20)
        data = dict(L=np.array([6,12,18,24]), sample_mean_Z=np.array([80.,81.,82.,83.]),
                    eta_C=np.array([.5,.4,.3,np.nan]), eta_C_sem=np.array([.1,.1,.1,0.]),
                    tau_Z=.02,tau_C=.3,accepted_L=np.nan,reference_L=24,n_images=20,
                    k_6=k,p_6=np.exp(-k),sem_6=np.zeros_like(k),
                    k_18=k,p_18=np.exp(-2*k),sem_18=np.zeros_like(k))
        np.savez_compressed(out/"custom.npz", **data)
        self.run_program("figures/regenerate_figure2_layout.py", "--image-file", image,
                         "--data-file", out/"custom.npz", "--output-dir", out/"fig2", "--pixel-size-mm", ".05", "--dpi", "35")
        audit = json.loads((out/"fig2/Figure_2_layout_audit.json").read_text())
        self.assertTrue(audit["identical_panel_boxes"])
        self.assertIsNone(audit["accepted_L"])
        self.assertEqual(audit["spectrum_supports"], [6,18])
        self.assertEqual(audit["n_images"], 20)
        self.assertEqual(audit["pixel_size_mm"], .05)

    def test_invalid_inputs(self):
        self.run_program("detect_stationary_domains_general.py", succeeds=False)
        self.run_program("property_full_homogenized_inside_property_masks.py",
                         "--kappa-maps", *self.kappa[:2], "--E-maps", *self.young[:1],
                         "--nu-maps", *self.poisson[:2], "--window-sizes", "4", "--L-rea", "4", succeeds=False)
        self.run_program("validation/rea_validation_wrapper.py", "--aggregate", "--checkpoint-dir",
                         self.work / "missing_checkpoints", succeeds=False)


if __name__ == "__main__":
    unittest.main()
