# setup.py
from setuptools import setup

setup(
    name='molstat',
    version='0.3',
    py_modules=['molstat'],
    install_requires=[
        'numpy',
        'matplotlib',
        'rdkit',
    ],
    entry_points={
        'console_scripts': [
            'molstat = molstat:main',
        ],
    },

    author='Zhaoqiang Chen',
    author_email='zqchen_simm@eqq.com',
    description='用于计算分子性质和绘制性质分布的工具。',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/autodataming/molstat',  # GitHub 仓库地址
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',



)

