CRPropa3
========

Development version of CRPropa. Use on your own risk.

Install from source
------------------------
1. Download the latest data archive from https://crpropa.desy.de/Interaction_data and extract to CRPropa3/data  
2. CRPropa uses CMAKE to configure. From the build directory call ccmake or cmake

        cd build
        ccmake ..  

3. Afterward configuring run make and make install as usual

        make
        make install 

The install path can be set with -DCMAKE_INSTALL_PREFIX=/my/path or with the option browser when using ccmake.  

Notes for Intel Compiler:  
use -DCMAKE_SHARED_LINKER_FLAGS="-lifcore" -DCMAKE_Fortran_COMPILER=ifort  

Provided Dependencies
---------------------
+ SOPHIA
    + for photo-hadronic interactions  
+ googletest 
    + for unit-tests  

Optional Dependencies
---------------------
+ Python and SWIG 
    + to use CRPropa from Python
    + tested for > Python 2.7
    + tested for > SWIG 2.0
+ FFTW3F 
    + for turbulent magnetic field grids
    + CRPropa needs the FFTW3 library compiled with the single precision option 
+ Gadget 
    + Magnetic fields for large scale structure data
+ OpenMP
    + for shared memory parallelization
+ googleperftools 
    + for performance optimizations regarding shared memory parallelization


XML Steering
------------
For backwards compatibility CRPropa 3 can be steered via XML cards ($cropra-xmlrun some_steeringcard.xml)
However, CRPropa 2 does not fully enforce the XML-1.0 standard (http://www.w3.org/XML).
To comply with the standard a few modifications to exisisting steering cards might have to be made.
Modification are  

1. XML-cards can have only one root node. The root node is

        <CRPropa>
        ...
        </CRPropa>

2. All nodes including the optional header node need to be closed, e.g.

        <?xml version="1.0" standalone=no?>
        <Option1> ... </Option1> or
        <Option2/>

3. Values need to be embraced by quotes, e.g. 

        <NumTrajectories="1000"/>