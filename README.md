# SmartFilter-API
API code for SmartFilter with Flasgger Docs.

## Requirements
* Docker  
* Docker Compose  

## Documentation
The `/apidocs/` page provides detailed information on every API call available.

* **Production:** [https://chiltepin.health.unm.edu/smartsfilter/apidocs/](https://chiltepin.health.unm.edu/smartsfilter/apidocs/)  
* **Local:** [http://localhost:8000/apidocs/](http://localhost:8000/apidocs/)  

**UI Repository:** [SmartFilter-UI](https://github.com/unmtransinfo/SmartFilter-UI/)

---

## Screenshots
Below are sample screenshots from the API UI:  

![API Home](/docs/images/api_docs.png)   

---

## Setup (Development)
1. `git clone https://github.com/yourusername/SmartFilter-API.git`  
2. `cd SmartFilter-API`  
3. Copy `.env.example` to `.env` (in the `/app` folder)  
4. Edit the `.env` credentials as needed  
5. Run:
   ```bash
   docker-compose -f compose-development.yml up --build


## Setup (Production on Chiltepin)
1. Copy `.env.example` to `.env`
2. Fill in/edit the `.env` credentials as needed 
3. Update apache2 config:
   - Create a new file for apache2 config: `/etc/apache2/sites-available/smartsfilterapi.conf`
   - Add the following line to `/etc/apache2/apache2.conf`:
     ```
     Include /etc/apache2/sites-available/smartsfilterapi.conf
     ```
   - Update the apache2 virtual config file: `/etc/apache2/sites-enabled/000-default.conf`
   - Run config check: `sudo apachectl configtest`
   - (If config check passed) reload apache: `sudo systemctl reload apache2`
4. (If server was previously up): `docker-compose -f compose-production.yml down`
5. Run `docker-compose -f compose-production.yml up --build -d`

## API Endpoints
1. /get_matchcounts
Description:
Counts matches of a single SMARTS pattern against a list of SMILES strings.

Query Parameters:

* smarts: The SMARTS string to match.

* SMILES: Comma or space-separated list of SMILES.

2. /get_matchFilter
Description:
Filters molecules that match a single SMARTS pattern.

Query Parameters:
* smarts: The SMARTS pattern to apply.
* SMILES: List of molecules to filter.

3. /get_multi_matchcounts
Description:
Counts matches of multiple SMARTS patterns against a list of SMILES strings.

Query Parameters:

* smarts: Comma/space-separated SMARTS patterns.

* SMILES: List of molecules to match.

4. /get_multi_matchfilter
Description:
Filters molecules that match all of the provided SMARTS patterns.

Query Parameters:
* smarts: Comma/space-separated SMARTS patterns.
* SMILES: List of molecules to filter.

5. /get_filterpains
Description:
Filters molecules against PAINS (Pan Assay Interference Compounds) substructure patterns.
Returns the filtered molecules and which PAINS patterns were matched.

Query Parameters:
* SMILES: List of SMILES strings to evaluate.

## Libraries Used

This API was developed with the help of the [rdkit-tools](https://github.com/jeremyjyang/rdkit-tools?tab=readme-ov-file#smarts) library for SMARTS pattern handling and molecular analysis.

**Thank you to the rdkit-tools developers** for providing this powerful library.

---

## Author / Support

**Author / Support:** Jeremy Yang, Jack Ringer  
**Developer:** Bivek Sharma Panthi

---

## Acknowledgment

This API was developed based on concepts and architecture inspired by the Badapple2-API repository:  
[https://github.com/unmtransinfo/Badapple2-API](https://github.com/unmtransinfo/Badapple2-API)

