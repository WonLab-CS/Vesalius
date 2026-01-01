{
  description = "Vesalius Development Shell";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-25.05";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        pkgs = import nixpkgs { inherit system; };

        rpkgs = builtins.attrValues {
          inherit (pkgs.rPackages)
              deldir
              sp
              tvR
              Matrix
              RColorBrewer
              future
              future_apply
              imagerExtra
              Signac
              ggpubr
              lmtest
              infix
              Seurat
              SeuratObject
              imager
              DESeq2
              RANN
              dplyr
              edgeR
              ggplot2
              igraph
              pwr
              purrr
              patchwork
              ggnewscale
              kohonen
              Rcpp
              RcppEigen
              TreeDist
              devtools
              tidyr;
        };

        
        
        system_packages = builtins.attrValues {
          inherit (pkgs)
            glibcLocales
            nix
            R
            pandoc
            which;
        };
      in
      {
        devShells.default = pkgs.mkShell {
        LOCALE_ARCHIVE = if pkgs.stdenv.hostPlatform.system == "x86_64-linux" then "${pkgs.glibcLocales}/lib/locale/locale-archive" else "";
        LANG = "en_US.UTF-8";
        LC_ALL = "en_US.UTF-8";
        LC_TIME = "en_US.UTF-8";
        LC_MONETARY = "en_US.UTF-8";
        LC_PAPER = "en_US.UTF-8";
        LC_MEASUREMENT = "en_US.UTF-8";

        buildInputs = rpkgs ++ system_packages;
        shellHook = ''
            echo "Vesalius Shell activated"

           '';
        };
      }
    );
}