/* ----------------------------------------------------------------------------
 * This file was automatically generated by SWIG (http://www.swig.org).
 * Version 3.0.12
 *
 * Do not make changes to this file unless you know what you are doing--modify
 * the SWIG interface file instead.
 * ----------------------------------------------------------------------------- */


public class fluidJNI {
 
  static { 
    try { 
        System.loadLibrary("fluid"); 
    } catch (UnsatisfiedLinkError e) { 
      System.err.println("Error while try to load the fluid shared library \n"); 
      System.exit(1); 
    } 
  } 

  public final static native void c_densitySolver(float[] jarg1, float[] jarg2, float jarg3, float[] jarg4, float[] jarg5, float jarg6, int jarg7, int jarg8);
  public final static native void c_velocitySolver(float[] jarg1, float[] jarg2, float[] jarg3, float[] jarg4, float[] jarg5, float[] jarg6, float jarg7, float jarg8, int jarg9, int jarg10);
}