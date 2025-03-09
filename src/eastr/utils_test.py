"""eastr utilities test."""

import unittest

from eastr import utils


class UtilsTest(unittest.TestCase):

    def test_sanitize_name(self):
        path = "/var/tmp/file.extension"
        result = utils.sanitize_name(path)
        expected = "file"
        self.assertEqual(result, expected)

    def test_sanitize_name_two(self):
        path = ".../file.extension"
        result = utils.sanitize_name(path)
        expected = "file"
        self.assertEqual(result, expected)

    def test_get_file_extension(self):
        path = ".../file.extension"
        result = utils.get_file_extension(path)
        expected = ".extension"
        self.assertEqual(result, expected)

    def test_sanitize_extension(self):
        path = ".../file.extension"
        extension = "_suffix.bam"
        result = utils.sanitize_and_update_extension(path, extension)
        expected = "file_suffix.bam"
        self.assertEqual(result, expected)


if __name__ == "__main__":
  unittest.main()
