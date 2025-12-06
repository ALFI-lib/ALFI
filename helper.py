#!/usr/bin/env python3
import argparse
import subprocess

dependencies = [
	'build-essential',
	'cmake',
	'git',
	'graphviz',
	'libgnuplot-iostream-dev',
	'libgtest-dev',
	'libqcustomplot-dev',
	'qtbase5-dev',
]

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--deps', action='store_true', help="Install dependencies")
parser.add_argument('-b', '--build', action='store_true', help="Build project (requires preset)")
parser.add_argument('-t', '--test', action='store_true', help='Run tests (requires preset)')
parser.add_argument('--doxygen', action='store_true', help='Generate Doxygen documentation')
parser.add_argument('preset', nargs='?', help="Preset to build and/or test (case-sensitive)")
args = parser.parse_args()

if (args.build or args.test) and not args.preset:
	parser.error("'--build' and '--test' require 'preset'")

if args.preset and not (args.build or args.test):
	parser.error("'preset' requires '--build' or '--test'")

def execute_command(command):
	print(' '.join(command), flush=True)
	subprocess.check_call(command)

if args.deps:
	execute_command(['sudo', 'apt', 'update'])
	execute_command(['sudo', 'apt', 'install', '-y'] + dependencies)
	execute_command(['git', 'submodule', 'update', '--init'])
if args.build:
	execute_command(['cmake', '--preset', args.preset])
	execute_command(['cmake', '--build', '--preset', args.preset, '-j'])
if args.test:
	execute_command(['ctest', '--preset', args.preset])
if args.doxygen:
	local = 'docs/doxygen/html/mathjax/es5/'
	remote = 'https://raw.githubusercontent.com/mathjax/MathJax/refs/tags/3.2.2/es5/'
	files = ['tex-svg.js', 'input/tex/extensions/physics.js']
	execute_command(['curl', '--create-dirs', '-C', '-', '-Z'] + [i for f in files for i in ['-o', local + f, remote + f]])
	execute_command(['doxygen', 'docs/doxygen/Doxyfile'])