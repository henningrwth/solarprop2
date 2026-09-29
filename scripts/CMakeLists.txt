set(SOLARPROP_SCRIPTS SOLARPROP
                      solarprop_arrows
                      solarprop_distributions
                      solarprop_trajectories)

install(PROGRAMS ${SOLARPROP_SCRIPTS} DESTINATION ${SKBUILD_SCRIPTS_DIR})
