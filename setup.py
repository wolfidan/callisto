from setuptools import setup, find_packages

def parse_requirements(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f if line.strip() and not line.startswith('#')]

setup(
    name='callisto',
    version='0.1.0',
    description='Set of scripts and functions to process data from the X-band CALLISTO spectrometer in Locarno Monti',
    author='Daniel Wolfensberger',
    author_email='daniel.wolfensberger@meteoswiss.com',
    packages=find_packages(include=['callisto_lib', 'callisto_lib.*']),
    python_requires='>=3.9',
    install_requires=parse_requirements('requirements.txt'),
    include_package_data=True,
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    entry_points={
        # Optional: create command-line tools from scripts
        # 'console_scripts': [
        #     'run-exp=scripts.run_experiment:main',
        # ],
    },
)