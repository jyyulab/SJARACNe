from setuptools import setup, find_packages

version = {}
with open("./version.py") as fp:
    exec(fp.read(), version)

setup(
    name="SJARACNe",
    version=version["__version__"],
    description="Gene network reverse engineering from big data",
    url="https://github.com/jyyulab/SJARACNe",
    author="Liang Ding, Chenxi Qian, Jiyang Yu",
    author_email="liang.ding@stjude.org, chenxi.qian@stjude.org, jiyang.yu@stjude.org",
    license="See LICENSE.md",
    install_requires=[
        "pandas >= 0.22.0",
        "numpy >= 1.14.2",
        "scipy == 1.0.1",
        "cwltool == 1.0",
    ],
    python_requires="==3.6.1",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    include_package_data=True,
    # test_suite="tests",
    entry_points={"console_scripts": ["sjaracne=SJARACNe.sjaracne:main"]},
)
