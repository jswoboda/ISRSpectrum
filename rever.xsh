$PROJECT = 'ISRSpectrum'
$ACTIVITIES = [
              'version_bump',  # Changes the version number in various source files (setup.py, __init__.py, etc)
              'changelog',  # Uses files in the news folder to create a changelog for release
              'tag',  # Creates a tag for the new version number
              'push_tag',  # Pushes the tag up to the $TAG_REMOTE
              'pypi',  # Sends the package to pypi
              #'conda_forge',  # Creates a PR into your package's feedstock
              'ghrelease'  # Creates a Github release entry for the new tag
               ]
$VERSION_BUMP_PATTERNS = [ ('setup.py', 'version\s*=.*,', "version='$VERSION',"),
('docs/source/conf.py', 'release\s*=.*,', "release=='$VERSION',", 'version\s*=.*,', "version='$VERSION',")]

$CHANGELOG_FILENAME = 'CHANGELOG.rst'  # Filename for the changelog
$CHANGELOG_TEMPLATE = 'TEMPLATE.rst'  # Filename for the news template

$PROJECT = "ISRSpectrum"
$WEBSITE_URL = "https://github.com/jswoboda/ISRSpectrum"
$GITHUB_ORG = "jswoboda"
$GITHUB_REPO = "ISRSpectrum"
$PYPI_NAME = "ISRSpectrum"
$PYPI_RC = '~/.pypirc'