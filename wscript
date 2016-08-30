#! /usr/bin/env python
# encoding: utf-8


def isxdigit(string):
	try:
		int(string, 16)
		return True
	except ValueError:
		return False

def check_output(*args):
	import subprocess
	process = subprocess.Popen(stdout = subprocess.PIPE, *args)
	output, unused_err = process.communicate()
	retcode = process.poll()
	return(output if retcode == 0 else None)

def get_hash(default_version_string):
	s = check_output(['git', 'rev-parse', 'HEAD'])
	return(s.decode().split('\n')[0] if s is not None else default_version_string)

def get_tag(hash):
	s = check_output(['git', 'show-ref', '--tags'])
	return(filter(lambda a: a == hash, s.decode().split('\n')) if s is not None else '')

def get_version_string(default_version_string):
	hash = get_hash(default_version_string)
	tag = get_tag(hash)
	hash = hash[0:7] if isxdigit(hash) == True else hash
	return('"%s"' % (tag if tag != '' else hash))


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
	# conf.env.append_value('CFLAGS', '-Wextra')
	# conf.env.append_value('CFLAGS', '-Wno-missing-field-initializers')
	# conf.env.append_value('CFLAGS', '-Wno-unused-parameter')
	conf.env.append_value('CFLAGS', '-Wno-unused-function')


	if conf.env.CC_NAME == 'icc':
		conf.env.append_value('CFLAGS', '-inline-max-size=20000')
		conf.env.append_value('CFLAGS', '-inline-max-total-size=50000')

	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-march=native')

	conf.env.append_value('LIBS', conf.env.LIB_Z + conf.env.LIB_BZ2 + conf.env.LIB_PTHREAD)
	conf.env.append_value('DEFINES', conf.env.DEFINES_Z + conf.env.DEFINES_BZ2 + ['COMB_VERSION_STRING=' + get_version_string("0.0.1")])
	conf.env.append_value('OBJS',
		['aw.o', 'fna.o', 'gaba_linear.o', 'gaba_affine.o', 'gaba_wrap.o', 'ggsea.o', 'gref.o', 'hmap.o', 'kopen.o', 'ngx_rbtree.o', 'psort.o', 'ptask.o', 'queue.o', 'queue_internal.o', 'sr.o', 'tree.o', 'zf.o'])


def build(bld):
	bld.recurse('arch')

	bld.objects(source = 'aw.c', target = 'aw.o')
	bld.objects(source = 'fna.c', target = 'fna.o')
	bld.objects(source = 'gaba.c', target = 'gaba_linear.o', defines = ['SUFFIX', 'MODEL=LINEAR'])
	bld.objects(source = 'gaba.c', target = 'gaba_affine.o', defines = ['SUFFIX', 'MODEL=AFFINE'])
	bld.objects(source = 'gaba_wrap.c', target = 'gaba_wrap.o')
	bld.objects(source = 'ggsea.c', target = 'ggsea.o')
	bld.objects(source = 'gref.c', target = 'gref.o')
	bld.objects(source = 'hmap.c', target = 'hmap.o')
	bld.objects(source = 'kopen.c', target = 'kopen.o')
	bld.objects(source = 'ngx_rbtree.c', target = 'ngx_rbtree.o')
	bld.objects(source = 'psort.c', target = 'psort.o')
	bld.objects(source = 'ptask.c', target = 'ptask.o')
	bld.objects(source = 'queue.c', target = 'queue.o')
	bld.objects(source = 'queue_internal.c', target = 'queue_internal.o')
	bld.objects(source = 'sr.c', target = 'sr.o')
	bld.objects(source = 'tree.c', target = 'tree.o')
	bld.objects(source = 'zf.c', target = 'zf.o')

	bld.program(
		source = ['comb.c'],
		target = 'comb',
		use = bld.env.OBJS,
		lib = bld.env.LIBS,
		defines = ['MAIN'],
		install_path = '${PREFIX}/bin')
