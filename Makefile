C_FLAGS := -O2

BUILD_DIR := build
DIST_DIR := dist
SRC_DIR := src

.PHONY: build
build: ${BUILD_DIR}/FastMap2 ${BUILD_DIR}/ExonWeaver ${BUILD_DIR}/MultiSeqNW

.PHONY: dist
dist: build
	mkdir -p dist/
	cp ${BUILD_DIR}/FastMap2 ${DIST_DIR}/
	cp ${BUILD_DIR}/ExonWeaver ${DIST_DIR}/
	cp ${BUILD_DIR}/MultiSeqNW ${DIST_DIR}/
	cp ${SRC_DIR}/*.pl ${DIST_DIR}/
	cp ${SRC_DIR}/*.pm ${DIST_DIR}/
	cp LICENSE ${DIST_DIR}/
	cp README.md ${DIST_DIR}/

.PHONY: check
check: build
	@echo "No checks yet"

.PHONY: clean
clean:
	rm -rf ${BUILD_DIR}
	rm -rf ${DIST_DIR}

${BUILD_DIR}/FastMap2: ${SRC_DIR}/FastMap2.c ${SRC_DIR}/BasicBio.c
	mkdir -p ${BUILD_DIR}
	${CC} ${C_FLAGS} -o $@ $^ -lm

${BUILD_DIR}/ExonWeaver: ${SRC_DIR}/ExonWeaver.c ${SRC_DIR}/BasicBio.c
	mkdir -p ${BUILD_DIR}
	${CC} ${C_FLAGS} -o $@ $^ -lm

${BUILD_DIR}/MultiSeqNW: ${SRC_DIR}/MultiSeqNW.c
	mkdir -p ${BUILD_DIR}
	${CC} ${C_FLAGS} -o $@ $^ -lm

# This Makefile knows how to build the project
# web site, which mostly consists of converting
# Markdown files to HTML.

PAGE_INPUTS := $(shell find docs -type f -name '*.md')
PAGE_OUTPUTS := $(PAGE_INPUTS:.md=.html)

HEADER := docs/media/header.html

SITE_TITLE := [Mirage]
PANDOC := pandoc -s -H ${HEADER}

.PHONY: docs
docs: ${PAGE_OUTPUTS}

.PHONY: serve
serve:
	@echo "Serving on http://localhost:8080"
	@echo "Use ctrl-c to terminate the server"
	@python3 -m http.server -d docs/ 8080

%.html: %.md ${HEADER}
	@echo "building $<"
	@${PANDOC} -M pagetitle:"/$(@D) ${SITE_TITLE}" -o "$@" "$<"
