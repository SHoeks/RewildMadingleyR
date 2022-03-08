\name{MadingleyRewilding-package}
\alias{MadingleyRewilding-package}
\alias{MadingleyRewilding}
\docType{package}

\title{MadingleyRewilding package}

\description{The MadingleyRewilding package is a modified version of the MadingleyR package intended to study the dynamics of rewilding actions.}

\details{
\strong{The MadingleyRewilding package has the following functions:}

\itemize{
  \item \code{MadingleyInit()}: Loads the default Madingley inputs
  \item \code{MadingleyInit()}: Initializes the Madingley model, includes some additional inputs parameters compared to the default run function of the MadingleyR package, see example link below
  \item \code{MadingleyRun()}: Runs the Madingley model, includes some additional inputs parameters compared to the default run function of the MadingleyR package, see example link below
  \item \code{DownloadLoadHalfDegreeInputs()}: Load the 0.5 degree Madingley spatial inputs, by default the MadingleyR package only comes with the 1 degree inputs
  \item \code{PlotCombinedTimelines()}: Plot the timelines of multiple consecutive simulations into a single plot
  \item \code{GetLibPath()}: Get the install folder of the package to copy the compiled C++ code when needed (when compiling the source code manually)
  \item \code{CropSpatialRastersToWindow()}: Crop spatial inputs to model domain, defined by spatial_window, helps to speed up loading times of consecutive MadingleyRewilding simulations

}

An example of how to setup and run a MadingleyRewilding simulation can be found on: \url{link here}\cr\cr


}
