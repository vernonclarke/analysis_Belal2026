import json
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

REPO_ROOT = Path(__file__).resolve().parents[1]
FUNCTIONS_DIR = REPO_ROOT / "Python functions"
sys.path.insert(0, str(FUNCTIONS_DIR))

from master_RNAscope import counts_from_roi_jsons, parse_field_metadata


class RNAscopeCountParsingTests(unittest.TestCase):
    def test_parse_field_metadata_accepts_dot_or_underscore_field_indices(self):
        self.assertEqual(
            parse_field_metadata("L1.ST8_C.UL_60x.01"),
            {
                "slice_id": "L1.ST8_C",
                "hemisphere": "UL",
                "field_index": 1,
                "condition": "Intact",
            },
        )
        self.assertEqual(
            parse_field_metadata("L4.ST8.B.L_60x07"),
            {
                "slice_id": "L4.ST8.B",
                "hemisphere": "L",
                "field_index": 7,
                "condition": "Lesioned",
            },
        )

    def test_parse_field_metadata_rejects_invalid_input(self):
        with self.assertRaises(ValueError):
            parse_field_metadata("")
        with self.assertRaises(ValueError):
            parse_field_metadata("not-a-field")

    def test_counts_from_roi_jsons_handles_none_and_empty_input(self):
        none_df = counts_from_roi_jsons(None)
        empty_df = counts_from_roi_jsons([])

        self.assertTrue(none_df.empty)
        self.assertTrue(empty_df.empty)
        self.assertEqual(list(none_df.columns), list(empty_df.columns))

    def test_counts_from_roi_jsons_rejects_invalid_field_name(self):
        with TemporaryDirectory() as tmpdir:
            json_path = Path(tmpdir) / "invalid.roi_analysis.json"
            json_path.write_text(
                json.dumps(
                    {
                        "field": "not-a-field",
                        "session": "L1.ST8",
                        "groups": {"NDNF+": [{"roi_index": 1, "dot_count": 3}]},
                    }
                )
            )

            with self.assertRaises(ValueError):
                counts_from_roi_jsons([json_path])

    def test_counts_from_roi_jsons_handles_large_input_batch(self):
        with TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            json_paths = []
            for index in range(200):
                field_index = (index % 20) + 1
                hemisphere = "UL" if index % 2 == 0 else "L"
                field = f"L1.ST8_C.{hemisphere}_60x.{field_index:02d}"
                json_path = root / f"{index:03d}.roi_analysis.json"
                json_path.write_text(
                    json.dumps(
                        {
                            "field": field,
                            "session": "L1.ST8",
                            "count_channel": "s_C004",
                            "groups": {
                                "NDNF+": [
                                    {"roi_index": 1, "dot_count": index},
                                    {"roi_index": 2, "dot_count": index + 1},
                                ]
                            },
                        }
                    )
                )
                json_paths.append(json_path)

            counts = counts_from_roi_jsons(json_paths)

        self.assertEqual(len(counts), 400)
        self.assertEqual(set(counts["condition"]), {"Intact", "Lesioned"})
        self.assertEqual(set(counts["count_channel"]), {"s_C004"})


if __name__ == "__main__":
    unittest.main()
