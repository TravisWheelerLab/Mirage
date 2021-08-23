C_FLAGS := -O2
PREFIX := /usr/local
WORK_DIR := dist

.PHONY: build
build: ${WORK_DIR}/FastMap2 ${WORK_DIR}/ExonWeaver ${WORK_DIR}/MultiSeqNW

.PHONY: install
install: build
	mkdir -p ${PREFIX}/bin
	cp ${WORK_DIR}/FastMap2 ${PREFIX}/bin/
	cp ${WORK_DIR}/ExonWeaver ${PREFIX}/bin/
	cp ${WORK_DIR}/MultiSeqNW ${PREFIX}/bin/

.PHONY: check
check: build
	@echo "No checks yet"

.PHONY: clean
clean:
	rm -rf ${WORK_DIR}

${WORK_DIR}/FastMap2: src/FastMap2.c src/BasicBio.c
	mkdir -p ${WORK_DIR}
	${CC} ${C_FLAGS} -o $@ $^ -lm

${WORK_DIR}/ExonWeaver: src/ExonWeaver.c src/BasicBio.c
	mkdir -p ${WORK_DIR}
	${CC} ${C_FLAGS} -o $@ $^ -lm

${WORK_DIR}/MultiSeqNW: src/MultiSeqNW.c
	mkdir -p ${WORK_DIR}
	${CC} ${C_FLAGS} -o $@ $^ -lm

