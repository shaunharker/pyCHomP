# installer script
rm -rf build
rm -rf dist
git submodule update --init --recursive
#pip install . --ignore-installed --no-cache-dir
pip install . --upgrade --no-deps --force-reinstall --no-cache-dir --user
