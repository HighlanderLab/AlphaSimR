#!/usr/bin/env Rscript

# Find the RcppTskit library directory and library file
# and return appropriate linker flags for use in Makevars
findRcppTskitLib <- function() {
  libdir <- system.file("libs", package = "RcppTskit")
  if (!nzchar(libdir)) {
    # nzchar() checks for non-zero character length
    stop("Unable to locate the RcppTskit shared library directory!")
  }

  libdirs <- libdir
  if (.Platform$OS.type == "windows") {
    r_arch <- sub("^/", "", .Platform$r_arch)
    if (nzchar(r_arch)) {
      libdirs <- c(libdirs, file.path(libdir, r_arch))
    } else {
      arch_dirs <- c("x64", "i386")
      arch_dirs <- arch_dirs[dir.exists(file.path(libdir, arch_dirs))]
      libdirs <- c(libdirs, file.path(libdir, arch_dirs))
    }
    libdirs <- unique(libdirs)
  }

  if (.Platform$OS.type == "unix") {
    candidates <- c(
      "RcppTskit.so", # Unix/Linux
      "RcppTskit.dylib" # macOS
    )
  } else if (.Platform$OS.type == "windows") {
    candidates <- c(
      "RcppTskit.dll.a", # MinGW/Rtools (default)
      "RcppTskit.lib", # MSVC (backup to MinGW/Rtools)
      "RcppTskit.dll" # DLL (backup to MinGW/Rtools)
    )
  } else {
    stop("Unknown .Platform$OS.type!")
  }
  libpaths <- unlist(
    lapply(libdirs, function(dir) file.path(dir, candidates)),
    use.names = FALSE
  )
  libfile <- libpaths[file.exists(libpaths)][1]
  if (is.na(libfile) || !nzchar(libfile)) {
    stop(
      "Unable to locate the RcppTskit library file in ",
      paste(libdirs, collapse = ", "),
      "!"
    )
  }

  if (.Platform$OS.type == "unix") {
    # Unix/Linux/macOS
    # Embed rpath on so the runtime loader can find RcppTskit.so.
    ret <- sprintf(
      "-Wl,-rpath,%s %s",
      shQuote(dirname(libfile)),
      shQuote(libfile)
    )
  } else if (.Platform$OS.type == "windows") {
    # Windows
    ret <- shQuote(libfile)
  } else {
    stop("Unknown .Platform$OS.type!")
  }
  return(ret)
}

# Render a Makevars file from a template by replacing placeholders.
# @param template character path to the template Makevars file
# @param output character path to the output Makevars file
renderMakevars <- function(template, output) {
  if (!file.exists(template)) {
    stop("Template file does not exist: ", template, "!")
  }
  RcppTskitLib <- findRcppTskitLib()
  lines <- readLines(con = template)
  lines <- gsub(
    x = lines,
    pattern = "@RCPPTSKIT_LIB@",
    replacement = RcppTskitLib,
    fixed = TRUE
  )
  writeLines(text = lines, con = output)
  invisible(TRUE)
}

# renderMakevars(template = "this_should_fail", output = "before_getting_to_output")
if (.Platform$OS.type == "unix") {
  # readLines(con = "src/Makevars.in")
  success <- renderMakevars(
    template = "src/Makevars.in",
    output = "src/Makevars"
  )
} else if (.Platform$OS.type == "windows") {
  # readLines(con = "src/Makevars.win.in")
  success <- renderMakevars(
    template = "src/Makevars.win.in",
    output = "src/Makevars.win"
  )
} else {
  stop("Unknown .Platform$OS.type!")
}
if (!success) {
  stop("renderMakevars() failed!")
}
