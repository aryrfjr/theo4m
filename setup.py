from setuptools import setup, find_packages
from pathlib import Path

this_dir = Path(__file__).parent
long_description = (this_dir / "README.md").read_text(encoding="utf-8")

setup(
    name="theo4m",
    version="0.1.0",
    description="Toolkit from Phys. Rev. Materials 4, 113603",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Ary Junior",
    packages=find_packages(),  # This will find the top-level theo4m/ package
    python_requires=">=3.7",
)

