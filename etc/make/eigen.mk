EIGEN_GIT_REPOSITORY ?= https://gitlab.com/libeigen/eigen.git
EIGEN_CMAKE = $(EIGEN_SRC_PATH)/CMakeLists.txt
EIGEN_MAKE = $(EIGEN_BUILD_PATH)/Makefile

EIGEN_MAIN = $(EIGEN_BUILD_PATH)/include/eigen3/Eigen/Eigen

$(EIGEN_CMAKE):
	mkdir -p $(@D)
	git clone -b $(EIGEN_BRANCH) $(EIGEN_GIT_REPOSITORY) $(@D)
	touch $@

${EIGEN_MAKE}: $(EIGEN_CMAKE)
	@echo ${EIGEN_MAIN} aus $(EIGEN_CMAKE)
	mkdir -p $(@D)
	cd $(@D) ; \
	cmake -DCMAKE_INSTALL_PREFIX=$(@D) $(<D)

${EIGEN_MAIN}: $(EIGEN_MAKE)
	cd ${<D} && $(MAKE) install
	test -f $(@) && touch $@

eigen: $(EIGEN_MAIN)
eigen-clean:
	rm -r $(EIGEN_BUILD_PATH)

EXTERNAL_DEPENDENCIES += eigen
INCLUDE_FLAGS += ${EIGEN_INCLUDE}
REQUIRED_MAKEVARS += EIGEN_BRANCH EIGEN_BUILD_PATH EIGEN_SRC_PATH EIGEN_GIT_REPOSITORY

.PHONY: eigen
