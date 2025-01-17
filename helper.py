#!/usr/bin/env python3
import argparse
import os
import subprocess
import sys

profile_directories = {
	'Debug': 'build/debug',
	'Release': 'build/release',
	'RelWithDebInfo': 'build/relwithdebinfo',
	'MinSizeRel': 'build/minsizerel',
	'Fast': 'build/fast',
	'Sanitize': 'build/sanitize',
}

dependencies = [
	'build-essential',
	'cmake',
	'libgnuplot-iostream-dev',
	'libgtest-dev',
	'libqcustomplot-dev',
	'qtbase5-dev',
]

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--deps', action='store_true', help="Install dependencies")
parser.add_argument('-b', '--build', action='store_true', help="Build project (requires profile)")
parser.add_argument('-t', '--test', action='store_true', help='Run tests (requires profile)')
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
		parser.error(f"Profile '{args.profile}' not found. Allowed profiles: '{"', '".join(profile_directories.keys())}'.")

try:
	if args.deps:
		apt_update_command = ['sudo', 'apt', 'update']
		print(' '.join(apt_update_command), flush=True)
		subprocess.check_call(apt_update_command)

		apt_install_command = ['sudo', 'apt', 'install', '-y'] + dependencies
		print(' '.join(apt_install_command), flush=True)
		subprocess.check_call(apt_install_command)

	if args.build:
		configure_command = ['cmake', '-DCMAKE_BUILD_TYPE=' + args.profile, '-B', profile_dir]
		print(' '.join(configure_command), flush=True)
		subprocess.check_call(configure_command)

		build_command = ['cmake', '--build', profile_dir, '-j', str(os.cpu_count())]
		print(' '.join(build_command), flush=True)
		subprocess.check_call(build_command)

	if args.test:
		test_command = ['ctest', '--test-dir', profile_dir, '--verbose']
		print(' '.join(test_command), flush=True)
		subprocess.check_call(test_command)
except subprocess.CalledProcessError as e:
	print(f'Error during build: {e}', file=sys.stderr)
	sys.exit(1)