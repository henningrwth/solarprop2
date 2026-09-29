set(MODEL_HEADERS
  burgerPotgieter.h
  jokipiiKopriva.h
  potgieterMoraal.h
  yamada.h
  standard2D.h
  helmod.h
  helmod2012.h
  helmod2018.h
  burger2012.h
  strauss2012.h
  custom.h
  )

set(MODEL_SOURCES
  burgerPotgieter.cc
  jokipiiKopriva.cc
  potgieterMoraal.cc
  yamada.cc
  standard2D.cc
  helmod.cc
  helmod2012.cc
  helmod2018.cc
  burger2012.cc
  strauss2012.cc
  custom.cc
  )

target_sources(solarprop PRIVATE ${MODEL_SOURCES} ${MODEL_HEADERS})

if(NOT DEFINED SKBUILD)
    target_sources(sprop_test PRIVATE ${MODEL_SOURCES} ${MODEL_HEADERS})
endif()
