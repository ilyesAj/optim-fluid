 /* fluid: SWIG configuration file to automatically generate the glue between C computing functions and java interface */
%include <arrays_java.i>
%apply float[] {float *};
%pragma(java) jniclasscode=%{ 
  static { 
    try { 
        System.loadLibrary("fluid"); 
    } catch (UnsatisfiedLinkError e) { 
      System.err.println("Error while try to load the fluid shared library \n"); 
      System.exit(1); 
    } 
  } 
%} 


 %module fluid
 %{
 /* Header to be exported in the JNI format by SWIG */
 extern void c_densitySolver(float *d, float *dOld, float diff, float *u, float *v , float dt, int grid_size, int vector_size);
 extern void c_velocitySolver( float *u, float *uOld, float *v, float *vOld, float *curl, float *d, float visc, float dt, int grid_size, int vector_size);
 %}
 
 extern void c_densitySolver(float *d, float *dOld, float diff, float *u, float *v , float dt, int grid_size, int vector_size);
 extern void c_velocitySolver( float *u, float *uOld, float *v, float *vOld, float *curl, float *d, float visc, float dt, int grid_size, int vector_size);









