"""eastr output utils test."""

import tempfile
import unittest

from eastr import output_utils


class OutputUtilsTest(unittest.TestCase):
  out_junctions_test_folder = None

  def setUp(self):
    self.out_junctions_test_folder = tempfile.TemporaryDirectory()

  def tearDown(self):
    self.out_junctions_test_folder.cleanup()

  def test_out_junctions_filelist_no_out_junctions(self):
    """Expect that function should return None if out_junctions in None."""
    results = output_utils.out_junctions_filelist(
      bam_list=None,
      gtf_path=None,
      bed_list=None,
      out_junctions=None,
      suffix="test"
    )
    self.assertIsNone(results)

  def test_out_junctions_filelist_gtf_path(self):
    """Expect that function should return a filelist for a gtf_path."""
    with tempfile.NamedTemporaryFile(suffix=".gtf") as gtfp:
      results = output_utils.out_junctions_filelist(
        bam_list=None,
        gtf_path=gtfp.name,
        bed_list=None,
        out_junctions=self.out_junctions_test_folder.name,
        suffix="test"
      )
      self.assertIsNotNone(results)
      self.assertTrue(results.endswith("test.bed"))

  def test_out_junctions_filelist_bedlist(self):
    """Expect that function should return a filelist for a bedlist."""
    # No results for suffix with "_original_junctions"
    results = output_utils.out_junctions_filelist(
      bam_list=None,
      gtf_path=None,
      bed_list=["test"],
      out_junctions=self.out_junctions_test_folder.name,
      suffix="_original_junctions"
    )
    self.assertIsNone(results)

   # No results for empty suffix"
    results = output_utils.out_junctions_filelist(
      bam_list=None,
      gtf_path=None,
      bed_list=["test"],
      out_junctions=self.out_junctions_test_folder.name,
      suffix=""
    )
    self.assertIsNone(results)

    # Results for valid bed_list and suffix
    sample_bed_list = ["test1", "test2"]
    results = output_utils.out_junctions_filelist(
      bam_list=None,
      gtf_path=None,
      bed_list=sample_bed_list,
      out_junctions=self.out_junctions_test_folder.name,
      suffix="_valid"
    )
    self.assertIsNotNone(results)
    for i, result in enumerate(results):
      self.assertTrue(result.endswith(f"{sample_bed_list[i]}_valid.bed"))

  def test_out_filtered_bam_filelist_none(self):
    """Expect that function should return a None for empty bamlist."""
    results = output_utils.out_filtered_bam_filelist(
      bam_list=None,
      out_filtered_bam="test",
      suffix=None,
    )
    self.assertIsNone(results)

  def test_out_filtered_bam_filelist_empty_out_filtered_bam(self):
    """Expect that function should return a None for empty out_filtered_bam."""
    results = output_utils.out_filtered_bam_filelist(
      bam_list=["test"],
      out_filtered_bam=None,
      suffix=None,
    )
    self.assertIsNone(results)

  def test_out_filtered_bam_filelist_single_bam_list(self):
    """Expect that function should return a single out filtered bam."""
    results = output_utils.out_filtered_bam_filelist(
      bam_list=["test"],
      out_filtered_bam=self.out_junctions_test_folder.name,
      suffix=None,
    )
    self.assertIsNotNone(results)
    self.assertEqual(len(results), 1)
    self.assertTrue(results[0].endswith("test_EASTR_filtered.bam"))

  def test_out_filtered_bam_filelist_single_bam_list_suffix(self):
    """Expect that function should return a single out a bam with suffix."""
    results = output_utils.out_filtered_bam_filelist(
      bam_list=["test"],
      out_filtered_bam=self.out_junctions_test_folder.name,
      suffix="_test",
    )
    self.assertIsNotNone(results)
    self.assertEqual(len(results), 1)
    self.assertTrue(results[0].endswith("test_test.bam"))

  def test_out_filtered_bam_filelist_multiple_bam_list_suffix(self):
    """Expect that function should return a multiple out a bam with suffix."""
    bam_list = ["test1", "test2"]
    results = output_utils.out_filtered_bam_filelist(
      bam_list=bam_list,
      out_filtered_bam=self.out_junctions_test_folder.name,
      suffix="_test",
    )
    self.assertIsNotNone(results)
    self.assertEqual(len(results), 2)
    for i, result in enumerate(results):
      self.assertTrue(result.endswith(f"{bam_list[i]}_test.bam"))

if __name__ == "__main__":
  unittest.main()
