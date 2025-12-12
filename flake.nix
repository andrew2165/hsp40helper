{
  description = "A Nix-flake-based RStudio Analysis environment";
  # use: nix develop

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-25.05";
  };

   outputs =
    { self, nixpkgs, ... }:
    let
      supportedSystems = [
        "aarch64-darwin"
        "x86_64-linux"
        "aarch64-linux"
      ];
      forAllSystems = nixpkgs.lib.genAttrs supportedSystems;
    in
    {
      devShells = forAllSystems (
        system:
        let
          pkgs = import nixpkgs { inherit system; };
        in
        {
          default = pkgs.mkShell {
            # create an environment with required R Packages
            packages = with pkgs; [
              (rWrapper.override { packages = [ rPackages.rstudio_prefs ]; })
              curlFull
              rPackages.curl
              rPackages.tidyverse
              rPackages.devtools
              rPackages.roxygen2
              rPackages.lintr
              rPackages.styler
              (rstudioWrapper.override {
                packages = with rPackages; [
                  tidyverse
                  drc
                  rstudio_prefs
		              devtools
                  roxygen2
                  lintr
                  styler
                ]; # add new R packages (from nix) here to get tied into Rstudio
              })
            ];

          shellHook = ''
            	      R -e "require(rstudio.prefs); rstudio_config_path('./rstudio-prefs.json');"
                    echo "RStudio Analysis Shell - ${system}"
            	      echo "Available commands: rstudio, R"
          '';
        };
    });

      # Also output rstudioPrefs for use in other flakes
      rstudioPrefs = ./rstudio-prefs.json;
    };
}
