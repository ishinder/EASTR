from __future__ import annotations as _annotations

import re
import subprocess

from mkdocs.config import Config
from mkdocs.structure.files import Files
from mkdocs.structure.pages import Page


def on_page_markdown(
    markdown: str, page: Page, config: Config, files: Files
) -> str:
  """Called on each file after it is read and before it is converted to HTML."""
  markdown = eastr_print_help(markdown, page)
  return markdown


def eastr_print_help(markdown: str, page: Page) -> str:
  if page.file.src_uri not in ("index.md"):
    return markdown

  output = subprocess.run(["eastr", "--help"], capture_output=True, check=True)
  logfire_help = output.stdout.decode()
  return re.sub(r"{{ *eastr_help *}}", logfire_help, markdown)
