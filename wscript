#! /usr/bin/env python
# encoding: utf-8


def isxdigit(string):
	try:
		int(string, 16)
		return True
	except ValueError, e:
		return False

def get_hash(default_version_string):
	import subprocess
	try:
		return(subprocess.check_output(['git', 'rev-parse', 'HEAD']).split('\n')[0])
	except subprocess.CalledProcessError, e:
		return(default_version_string)

def get_tag(hash):
	import subprocess
	try:
		return(filter(lambda a: a == hash, subprocess.check_output(['git', 'show-ref', '--tags']).split('\n')))
	except subprocess.CalledProcessError, e:
		return('')

def get_version_string(default_version_string):
	hash = get_hash(default_version_string)
	tag = get_tag(hash)
	hash = hash[0:7] if isxdigit(hash) == True else hash
	return('"{}"'.format(tag if tag != '' else hash))


def options(opt):
	opt.load('compiler_c')
	opt.recurse('arch')

def configure(conf):
	conf.load('ar')
	conf.load('compiler_c')

	conf.recurse('arch')

	if 'LIB_Z' not in conf.env:
		conf.check_cc(
			lib = 'z',
			defines = ['HAVE_Z'],
			mandatory = False)

	if 'LIB_BZ2' not in conf.env:
		conf.check_cc(
			lib = 'bz2',
			defines = ['HAVE_BZ2'],
			mandatory = False)

	if 'LIB_PTHREAD' not in conf.env:
		conf.check_cc(lib = 'pthread')

	conf.env.append_value('CFLAGS', '-Wall')
	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-Wno-unused-function')
	conf.env.append_value('CFLAGS', '-march=native')
	conf.env.append_value('CFLAGS', '-flto')
	conf.env.append_value('LDFLAGS', '-flto')

	conf.env.append_value('LIB_COMB', conf.env.LIB_Z + conf.env.LIB_BZ2 + conf.env.LIB_PTHREAD)
	conf.env.append_value('DEFINES_COMB', conf.env.DEFINES_Z + conf.env.DEFINES_BZ2 + ['COMB_VERSION_STRING=' + get_version_string("0.0.1")])
	conf.env.append_value('OBJ_COMB',
		['aw.o', 'fna.o', 'gaba_linear.o', 'gaba_affine.o', 'gaba_wrap.o', 'ggsea.o', 'gref.o', 'hmap.o', 'kopen.o', 'ngx_rbtree.o', 'psort.o', 'ptask.o', 'queue.o', 'queue_internal.o', 'sr.o', 'tree.o', 'zf.o'])


def build(bld):
	bld.recurse('arch')

	bld.objects(source = 'aw.c', target = 'aw.o')
	bld.objects(source = 'fna.c', target = 'fna.o')
	bld.objects(source = 'gaba.c', target = 'gaba_linear.o', includes = ['.'], defines = ['SUFFIX', 'MODEL=LINEAR'])
	bld.objects(source = 'gaba.c', target = 'gaba_affine.o', includes = ['.'], defines = ['SUFFIX', 'MODEL=AFFINE'])
	bld.objects(source = 'gaba_wrap.c', target = 'gaba_wrap.o', includes = ['.'])
	bld.objects(source = 'ggsea.c', target = 'ggsea.o')
	bld.objects(source = 'gref.c', target = 'gref.o', includes = ['.'])
	bld.objects(source = 'hmap.c', target = 'hmap.o')
	bld.objects(source = 'kopen.c', target = 'kopen.o')
	bld.objects(source = 'ngx_rbtree.c', target = 'ngx_rbtree.o')
	bld.objects(source = 'psort.c', target = 'psort.o', includes = ['.'])
	bld.objects(source = 'ptask.c', target = 'ptask.o')
	bld.objects(source = 'queue.c', target = 'queue.o')
	bld.objects(source = 'queue_internal.c', target = 'queue_internal.o')
	bld.objects(source = 'sr.c', target = 'sr.o')
	bld.objects(source = 'tree.c', target = 'tree.o')
	bld.objects(source = 'zf.c', target = 'zf.o')

	bld.program(
		source = ['comb.c'],
		target = 'comb',
		use = bld.env.OBJ_COMB,
		lib = bld.env.LIB_COMB,
		defines = ['MAIN'] + bld.env.DEFINES_COMB)
