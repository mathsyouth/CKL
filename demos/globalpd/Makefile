# Must set these gl and glut locations to build 'falling'

GL_INCPATH = -I/usr/local/include/
GL_LIBPATH = -L/usr/local/lib/ -L/usr/X11R6/lib/
GL_LIBS  =  -lglut -lGLU -lGL -lXext -lXmu -lXi -lX11

XML_INCPATH = -I./tinyxml
XML_LIBPATH = -L./tinyxml
XML_LIBS = -ltinyxml


.SUFFIXES: .cpp

CC = g++
CFLAGS  = -O2 -I. -I../../include $(GL_INCPATH) $(XML_INCPATH)
LDFLAGS	= -L. -L../../lib $(GL_LIBPATH) $(XML_LIBPATH)
LDLIBS  = -lCKL -lm       $(GL_LIBS) $(XML_LIBS)

SRCS    = main.cpp 

OBJECTS	= main.o

TARGET = globalpd

CLEAN   = $(OBJECTS) $(TARGET)

.cpp.o:
	$(CC) ${CFLAGS} -c $<

$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJECTS) -L. $(LDFLAGS) $(LDLIBS)

run: $(TARGET)
	$(TARGET)

clean:
	/bin/rm -f $(CLEAN)

