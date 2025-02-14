#!/usr/bin/env python3
import argparse
import os
import subprocess

profile_directories = {
	'Debug': 'build/debug',
	'Release': 'build/release',
	'RelWithDebInfo': 'build/relwithdebinfo',
	'MinSizeRel': 'build/minsizerel',
	'Fast': 'build/fast',
	'FastParallel': 'build/fastparallel',
	'Sanitize': 'build/sanitize',
}

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
parser.add_argument('-b', '--build', action='store_true', help="Build project (requires profile)")
parser.add_argument('-t', '--test', action='store_true', help='Run tests (requires profile)')
parser.add_argument('--doxygen', action='store_true', help='Generate Doxygen documentation')
parser.add_argument('profile', nargs='?', help="Profile to build and/or test (case-insensitive)")
args = parser.parse_args()

if (args.build or args.test) and not args.profile:
	parser.error("'--build' and '--test' require 'profile'")

if args.profile and not (args.build or args.test):
	parser.error("'profile' requires '--build' or '--test'")

profile_dir = None
if args.profile:
	try:
		profile_dir = next(d for p, d in profile_directories.items() if p.lower() == args.profile.lower())
	except StopIteration:
		parser.error(f"Profile '{args.profile}' not found. Allowed profiles: '{"', '".join(profile_directories)}'.")

def execute_command(command):
	print(' '.join(command), flush=True)
	subprocess.check_call(command)

if args.deps:
	execute_command(['sudo', 'apt', 'update'])
	execute_command(['sudo', 'apt', 'install', '-y'] + dependencies)
	execute_command(['git', 'submodule', 'update', '--init'])
if args.build:
	execute_command(['cmake', '-DCMAKE_BUILD_TYPE=' + args.profile, '-B', profile_dir])
	execute_command(['cmake', '--build', profile_dir, '-j', str(os.cpu_count())])
if args.test:
	execute_command(['ctest', '--test-dir', profile_dir, '--verbose'])
if args.doxygen:
	local = 'docs/doxygen/html/mathjax/es5/'
	remote = 'https://raw.githubusercontent.com/mathjax/MathJax/refs/tags/3.2.2/es5/'
	files = ['tex-svg.js', 'input/tex/extensions/physics.js']
	execute_command(['curl', '--create-dirs', '-C', '-', '-Z'] + [i for f in files for i in ['-o', local + f, remote + f]])
	execute_command(['doxygen', 'docs/doxygen/Doxyfile'])