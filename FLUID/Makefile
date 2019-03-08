
CC=gcc
INCLUDE=-I/usr/lib/jvm/java-8-openjdk-amd64/include/ -I/usr/lib/jvm/java-8-openjdk-amd64/include/linux/
CFLAGS=-c -fpic $(INCLUDE)
IFLAGS=-g fluid.c -qopt-report-annotate=html 
BLASFLAGS=-I /opt/OpenBLAS/include/ 
BLASLIB= -lopenblas -lpthread -lgfortran
OPENMPFLAGS=-fopenmp
all: core java
optim1: coreopt1 java
optim2: coreopt2 java
optim3: coreopt3 java

coreopt3:  interface_c_java_wrap.o fluidopt3.o
	$(CC) $(OPENMPFLAGS) -shared interface_c_java_wrap.o fluidopt3.o -o libfluid.so
coreopt2:  interface_c_java_wrap.o fluidopt1.o
	$(CC) -ffast-math -shared interface_c_java_wrap.o fluidopt1.o -o libfluid.so
coreopt1:  interface_c_java_wrap.o fluidopt1.o
	$(CC) -O3 -shared interface_c_java_wrap.o fluidopt1.o -o libfluid.so
core:  interface_c_java_wrap.o fluid.o
	$(CC) -shared interface_c_java_wrap.o fluid.o -L/opt/OpenBLAS/lib $(BLASLIB) $(OPENMPFLAGS) -o libfluid.so

fluid.o: fluid.c
	$(CC) $(CFLAGS) $(BLASLIB) $(BLASFLAGS) $(OPENMPFLAGS) $<

interface_c_java_wrap.o: 
	swig -java interface_c_java.swig
	$(CC) $(CFLAGS) $(BLASLIB) $(BLASFLAGS) $(OPENMPFLAGS) interface_c_java_wrap.c
java:
	javac *.java

run:
	appletviewer  -J"-Djava.security.policy=applet.policy" demo.html 

debug:
	appletviewer -debug -J"-Djava.security.policy=applet.policy" demo.html 

fluidopt1.o: fluidopt1.c
	$(CC) $(CFLAGS) $<
fluidopt3.o: fluidopt3.c
	$(CC) $(CFLAGS) $(OPENMPFLAGS) $<
clean:
	rm -f *.o
	rm -f *.so
	rm -f fluid.java
	rm -f fluidJNI.java
	rm -f SWIGTYPE_p_float.java
	rm -f interface_c_java_wrap.c
	rm -f *.class
