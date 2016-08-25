
all: build

configure:
	python waf configure

build: configure
	python waf build

clean: configure
	python waf clean
