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

