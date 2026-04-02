# HyPlAs C++ - Makefile
# ======================
#
# Build: ./configure && make
# Install: make install
#

# ============================================================================
# Configuration
# ============================================================================

# Include generated config if available, otherwise use defaults
-include config.mk

# Default values (used if configure hasn't been run)
CXX      ?= c++
CXXFLAGS ?= -std=c++20 -Wall -Wextra -Wpedantic -O2 -DNDEBUG
LDFLAGS  ?= -lz
PREFIX   ?= /usr/local
BINDIR   ?= $(PREFIX)/bin

# Directories
SRCDIR   := src
BUILDDIR := build

# ============================================================================
# hyplas - unified binary with pipeline and utilities
# ============================================================================
HYPLAS_TARGET := $(BUILDDIR)/hyplas
HYPLAS_SRCS := $(SRCDIR)/hyplas_main.cpp \
               $(SRCDIR)/process.cpp \
               $(SRCDIR)/pipeline.cpp
HYPLAS_OBJS := $(patsubst $(SRCDIR)/%.cpp,$(BUILDDIR)/%.o,$(HYPLAS_SRCS))

# All sources for dependency tracking
ALL_SRCS := $(HYPLAS_SRCS)
ALL_OBJS := $(HYPLAS_OBJS)

# Headers
HEADERS := $(wildcard $(SRCDIR)/*.hpp $(SRCDIR)/*.h)

# ============================================================================
# Targets
# ============================================================================

.PHONY: all
all: $(HYPLAS_TARGET)

# Create build directory
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# hyplas unified binary
$(HYPLAS_TARGET): $(HYPLAS_OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(HYPLAS_OBJS) $(LDFLAGS)

# Compile any .cpp to .o in build directory
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(HEADERS) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# ============================================================================
# Install
# ============================================================================

.PHONY: install
install: all
	@echo "Installing to $(BINDIR)"
	@mkdir -p "$(BINDIR)"
	cp "$(HYPLAS_TARGET)" "$(BINDIR)/"
	@echo "Installation complete."

.PHONY: uninstall
uninstall:
	@echo "Uninstalling from $(BINDIR)"
	rm -f "$(BINDIR)/hyplas"
	@echo "Uninstall complete."

# ============================================================================
# Clean
# ============================================================================

.PHONY: clean
clean:
	rm -rf $(BUILDDIR)

.PHONY: distclean
distclean: clean
	rm -f config.mk

.PHONY: rebuild
rebuild: clean all

# ============================================================================
# Development utilities
# ============================================================================

# Generate compile_commands.json for LSP/clangd
.PHONY: compile_commands
compile_commands:
	@echo '[' > compile_commands.json
	@first=1; for src in $(ALL_SRCS); do \
		[ $$first -eq 0 ] && echo ',' >> compile_commands.json; \
		first=0; \
		obj=$(BUILDDIR)/$$(basename $${src%.cpp}.o); \
		echo '{"directory":"$(CURDIR)","command":"$(CXX) $(CXXFLAGS) -c -o '$$obj' '$$src'","file":"'$$src'"}' >> compile_commands.json; \
	done
	@echo ']' >> compile_commands.json
	@echo "Generated compile_commands.json"

# ============================================================================
# Info / Help
# ============================================================================

.PHONY: info
info:
	@echo "Build configuration:"
	@echo "  CXX:           $(CXX)"
	@echo "  CXXFLAGS:      $(CXXFLAGS)"
	@echo "  LDFLAGS:       $(LDFLAGS)"
	@echo ""
	@echo "Installation paths:"
	@echo "  PREFIX:        $(PREFIX)"
	@echo "  BINDIR:        $(BINDIR)"
	@echo ""
	@echo "Target:"
	@echo "  HYPLAS_TARGET: $(HYPLAS_TARGET)"

.PHONY: help
help:
	@echo "HyPlAs C++ build system"
	@echo ""
	@echo "Quick start:"
	@echo "  ./configure            Configure with defaults (auto-detects env)"
	@echo "  ./configure --prefix=PATH  Configure with custom prefix"
	@echo "  make                   Build hyplas binary"
	@echo "  make install           Install to configured prefix"
	@echo ""
	@echo "Commands:"
	@echo "  hyplas run [options]              Run assembly pipeline"
	@echo "  hyplas innotin [options]          Filter reads not in subsets"
	@echo "  hyplas select-missing-reads ...   Extract reads from PAF"
	@echo "  hyplas split-plasmid-reads ...    Classify reads"
	@echo ""
	@echo "Build targets:"
	@echo "  all              Build hyplas binary (default)"
	@echo ""
	@echo "Installation:"
	@echo "  install          Install binary to BINDIR"
	@echo "  uninstall        Remove installed binary"
	@echo ""
	@echo "Maintenance:"
	@echo "  clean            Remove build artifacts"
	@echo "  distclean        Remove build artifacts and config.mk"
	@echo "  rebuild          Clean and rebuild all"
	@echo ""
	@echo "Development:"
	@echo "  compile_commands Generate compile_commands.json for LSP"
	@echo "  info             Print build configuration"
	@echo "  help             Show this message"
	@echo ""
	@echo "Run './configure --help' for configuration options."
