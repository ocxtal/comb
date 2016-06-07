#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.recurse('sr')
	opt.recurse('aw')
	opt.recurse('gref')
	opt.recurse('ggsea')
	opt.recurse('ptask')
	opt.load('compiler_c')

def configure(conf):
	conf.recurse('sr')
	conf.recurse('aw')
	conf.recurse('gref')
	conf.recurse('ggsea')
	conf.recurse('ptask')

	conf.load('ar')
	conf.load('compiler_c')

	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-std=c99')
	conf.env.append_value('CFLAGS', '-march=native')

	conf.env.append_value('LIB_COMB',
		conf.env.LIB_SR + conf.env.LIB_AW + conf.env.LIB_GREF + conf.env.LIB_GGSEA + conf.env.LIB_PTASK)
	conf.env.append_value('DEFINES_COMB',
		conf.env.DEFINES_SR + conf.env.DEFINES_AW + conf.env.DEFINES_GREF + conf.env.DEFINES_GGSEA + conf.env.DEFINES_PTASK)
	conf.env.append_value('OBJ_COMB',
		conf.env.OBJ_SR + conf.env.OBJ_AW + conf.env.OBJ_GREF + conf.env.OBJ_GGSEA + conf.env.OBJ_PTASK)


def build(bld):
	bld.recurse('sr')
	bld.recurse('aw')
	bld.recurse('gref')
	bld.recurse('ggsea')
	bld.recurse('ptask')

	bld.program(
		source = ['comb.c'],
		target = 'comb',
		use = bld.env.OBJ_COMB,
		lib = bld.env.LIB_COMB,
		defines = ['MAIN'] + bld.env.DEFINES_COMB)
