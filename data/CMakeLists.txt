set(SOLARPROP_DATA angle.dat
                   nm.dat
                   omniweb.dat
                   ssn.dat
		   )

install(FILES ${SOLARPROP_DATA} DESTINATION ${SKBUILD_PLATLIB_DIR}/solarprop-data)
