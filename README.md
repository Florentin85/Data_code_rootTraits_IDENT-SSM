# Explanation of variable names for root trait data from the International Diversity Experiment Network site in Sault Ste. Marie (Ontario, Canada, IDENT-SSM)
## Usage
The code and data provided in this repository are intended for research and educational purposes. If you use any part of this repository in your work, please cite the original research paper:
```bibtex
@article{jaegerYoungTemperateTree2025,
  title = {Young Temperate Tree Species Show Fine Root Trait Acclimation to Differences in Water Availability},
  author = {Jaeger, Florentin C. and Messier, Christian and Beyer, Friderike and Aubin, Isabelle and Parker, William C. and Handa, I. Tanya},
  date = {2025-11-10},
  journaltitle = {Plant and Soil},
  shortjournal = {Plant Soil},
  issn = {1573-5036},
  doi = {10.1007/s11104-025-08047-5},
  url = {https://doi.org/10.1007/s11104-025-08047-5},
  urldate = {2025-11-11},
  abstract = {Growing-season water availability in forests is being altered due to climate change. We examined how the magnitude and plasticity of absorptive fine root traits of six young temperate tree species responded to high versus low soil water availability over five years.},
  langid = {english},
  keywords = {drought,functional traits,morphology,plasticity,root architecture,root traits},
}
```

### How to Download and Clone the Repository

1. **Clone the Repository**:
   To clone this repository, you will need Git installed on your local machine. Open your terminal or command prompt and run the following command:

   ```bash
   git clone https://github.com/Florentin85/Data_code_rootTraits_IDENT-SSM.git
   ```

2. **Run the Analysis**:
  You can open the R Quarto files and run the analysis directly. To run a specific R Quarto file, use the following R command in your R environment:

    ```r
    quarto::quarto_render("your_file.qmd")
    ```
  Replace your_file.qmd with the actual name of the R Quarto file you want to render.

3. **Explore the Data**:
  The data files are included in the directory in .csv format. You can load and explore them using standard R functions, such as read.csv() or readRDS().

4. **Reproduce the Analysis**:
  Follow the instructions in the R Quarto files to reproduce the analysis. Each file contains detailed comments and explanations of the code and steps.

5. **Report Issues**:
  If you encounter any issues or have questions, please open an issue on the GitHub repository or contact me directly at the email provided below.

# Methodology and variables explained
| **Explanation** | **Name** |
| ----------- | ---- |
Identification rock volume | ID_rock
tree species | sp
Resolution in WinRhizo set to 500 dpi | resolution
day of year 2018 | day of year
date of field sampling | sampling date
root fraction (absorptive fine roots AF, transportive fine roots TF, coarse roots C, dead roots D, fragments F (<1 cm length) | fraction
fresh weight gram with 100 % fragments added to AF and D (F multiplied by 10 to get from 10 % to 100 % root fragments and distributed proportionally according to the fractions sorted by species and live and dead fine roots.) | fresh weight aff
partial fresh weight gram | part fresh w
dry weight gram | dry weight
dry weight gram fragments 100 % (F x 10) | dry weight f 100
calculated total dry weight in grams (g) | calc tot dry weight
total dry weight gram with F = 100 %, same method as for fresh weight aff. Excel formula: f100 / (AF + D) x AF + AF = af proportional (same formula according to D --> … (D + … ) x D + D | tot dry weight aff
Percentage fragments on tot dry weight aff | f perc prop weight aff
volume calculated rocks in milliliters (ml) per sample (strata: 5 cm and 10 cm) | vol calc sk ml
% skeleton in soil corer per stratum --> (vol calc sk ml / soil vol cm^-3^) x 100 | vol perc sk
vol fine soil cm^-3^ -->  volume corer - volume skeleton | vol fine cm^-3^
day of dry roots weighing | date dry weighing
Soil volume in m^-3^ from the root corer, diameter 5.08 cm and height of soil cores 5 and 10 cm respectively (0-5, 5-10, 10-15, 15-20 & 20-30 cm depth) | soil vol m^-3^
Fine root biomass = tot dry weight aff in gram / volume fine soil in cubic centimeter (Fine root biomass is fine root density in this case, biomass is related to m^-2^ and not cm^-3^) | FRB g cm^-3^
Specific root length = meters root length aff 100 cm divided by 100 (to get meters) divided by total dry weight aff (gram) | SRL m g^-1^
Root length density = (length aff cm / 100) / (vol fine cm^-3^ / 1000000) One million to get m^-3^. | RLD m m^-3^
We noted, if there was mycelium present in the scan --> Yes/No (Changed grey scale to get best scan result!) | mycelium
length from WinRhizo (total length per fraction, each fraction was scanned individually) | length cm
length from WinRhizo for fractions C & TF. F multiplied by 10 to get 100 % fine root fragments. Formula: f100 / (AF + D) x AF + AF --> to get length from F100 % proportionally distributed according to AF & D. Same for D. | length aff cm
Average root volume from WinRhizo --> Wrong calculation!!! --> WinRhizo assumes roots are perfect cylinders! | rootvol cm^-3^
Sum of the root volume per diameter class from WinRhizo in cubic centimeter (10 classes, 5 mm each) | rootvol class cm^-3^
Sum root volume for diameter classes. Volume from WinRhizo for fractions C & TF. F multiplied by 10 to get 100 % fine root fragments. Formula: f100 / (AF + D) x AF + AF --> to get volume from F100 % proportionally distributed according to AF & D. Same for D. | rootvol aff class cm^-3^
Root Tissue Density (RTD) --> Total dry weight aff / root vol aff class cm^-3^ | rtd g cm^-3^
number of root tips per root length (length aff cm). Fragments **NOT** multiplied by 10 and **NOT** distributed proportionally according to AF and D, because WinRhizo would consider fragments ends as root tips and this would be unrealistic! (same with forks) | ntips n length cm^-1^
Root dry matter content (RDMC) dry g / fresh g -->  total dry weight aff / fresh weight aff | rdmc dry g^-1^ fresh g^-1^

# License
This repository is licensed under the MIT License - see the LICENSE file for details.

# Contact
For any questions or feedback, please contact: jaeger.florentin_clemens[at]courier.uqam.ca
