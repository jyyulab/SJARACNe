#############################################################################
# Makefile for building: sjaracne
# Template: app
#############################################################################

UNAME_S := $(shell uname -s)

# CC
ifeq ($(UNAME_S), Darwin)
	CC := clang++ -arch x86_64
else
	CC := g++
endif

# Folders
SRCDIR := src
BUILDDIR := build
TARGETDIR := bin

# Targets
EXECUTABLE := sjaracne.exe
TARGET := $(TARGETDIR)/$(EXECUTABLE)

# Final Paths
INSTALLBINDIR := /usr/local/bin

# Code Listss
SRCEXT := cpp
HEADEXT := h
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%, $(BUILDDIR)/%, $(SOURCES:.$(SRCEXT)=.o))

# Folder Lists
INCDIRS := $(shell find src/* -name '*.h' -exec dirname {} \; | sort | uniq)
INCLIST := $(patsubst src/%, -I src/%, $(INCDIRS))
BUILDLIST := $(patsubst src/%, $(BUILDDIR)/%, $(INCDIRS))

# Shared Compiler Flags
CFLAGS := -c -pipe -g -Wall -W -fPIC
INC := -I src $(INCLIST) -I /usr/local/include
LIB := -L /usr/local/lib

ifeq ($(UNAME_S), Linux)
	CFLAGS += -O3
else
	CFLAGS += -std=c++11 -stdlib=libc++ -O3
endif

$(TARGET): $(OBJECTS)
	@mkdir -p $(TARGETDIR)
	@echo "Linking..."
	@echo "		Linking $(TARGET)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDLIST) $(BUILDDIR)
	@echo "Compiling $<..."; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

clean:
	@echo "Cleaning $(TARGET)..."; $(RM) -r $(BUILDDIR)/* $(TARGET)

install:
	@echo "Installing $(EXECUTABLE)..."; cp $(TARGET) $(INSTALLBINDIR)

distclean:
	@echo "Removing $(EXECUTABLE)"; rm $(INSTALLBINDIR)/$(EXECUTABLE)
