# installer script
rm -rf build
git submodule update --init --recursive
#pip install . --ignore-installed --no-cache-dir
pip install . --upgrade --no-deps --force-reinstall --no-cache-dir
