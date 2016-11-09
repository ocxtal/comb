
CC=gcc
PREFIX=/usr/local

all: build

configure:
	python waf configure CC=${CC} --prefix=${PREFIX}

build: configure
	python waf build

clean: configure
	python waf clean

install:
	python waf configure CC=${CC} --prefix=${PREFIX} install

