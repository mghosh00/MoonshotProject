import sqlite3
import re
import csv
from pathlib import Path


class DatabaseManager:
    def __init__(self, database_path: str):
        """
        Parameters
        ----------
        database_path: str
            Full path of existing SQLite database or one to be created automatically
        """

        # Store SQLite database path for reference
        self._database_path = database_path

        # Connect to and Create (if doesn't already exist) SQLite database
        self._conn = sqlite3.connect(database_path)

        # List of tables that need to be dropped to reset the SQLite database
        self._drop_order = ('compounds', 'submissions', 'assays',)

    def get_conn(self):
        """
        get_conn returns a connection to the SQLite database self._database_path

        Returns
        -------
            SQLite connection object
        """
        return self._conn

    def drop_all(self):
        """
        drop_all drops all tables created by this class to reset the SQLite database
        """

        # Get connection to SQLite database
        conn = self.get_conn()

        # Drop all tables in dependency order
        for table_name in self._drop_order:
            conn.execute('DROP TABLE IF EXISTS ' + table_name)

    def create(self):
        """
        create - creates all tables required by this class in the SQLite database
        """

        # Get connection to SQLite database
        conn = self.get_conn()

        # Create a table to store all COVID Moonshot submissions
        conn.execute('''
CREATE TABLE submissions
(
    submission_id VARCHAR(20) PRIMARY KEY,
    name_code CHAR(3) not null,
    institute_code CHAR(3) not null,
    random_id CHAR(8) not null
)
        ''')

        # Create a table to store all COVID Moonshot compound submissions
        conn.execute('''
CREATE TABLE compounds
(
    compound_id VARCHAR(20) PRIMARY KEY,
    smiles VARCHAR(2000) not null,
    submission_id VARCHAR(20) not null,
    made VARCHAR(5) not null,
    inchi_key CHAR(27) not null,
    MW DECIMAL not null,
    r_avg_IC50 DECIMAL,
    FOREIGN KEY(submission_id) REFERENCES submission_id(submission_id) 
)
        ''')

        # Create a table to store specific assays
        conn.execute('''
CREATE TABLE assays
(
    id UNIQUEIDENTIFIER PRIMARY KEY,
    compound_id VARCHAR(20) not null,
    f_avg_IC50 DECIMAL,
    r_avg_IC50 DECIMAL,
    covalent_fragment boolean,
    covalent_warhead boolean,
    FOREIGN KEY(compound_id) REFERENCES compounds(compound_id)
)
        ''')

    def populate_submissions_table(self, all_data_file: Path):
        """Populate the table submissions by reading out all of the unique submissions from $all_data_file
        """

        # Get database connection
        conn = self.get_conn()

        # Compile REGEX for submission ID
        sub_pattern = re.compile('(\w{3})-(\w{3})-(\w{8})-(\d+)')

        # Open the data file
        with open(all_data_file, 'r', encoding="utf8") as fh:
            # Initialise a CSV reader
            reader = csv.reader(fh, delimiter=',')

            # Create dictionary to generate unique list with
            unique_submissions = {}

            # Iterate rows in data file
            for cols in reader:
                # Get submission ID
                sub_id = cols[1]

                # Apply regular expression
                # Lots of ways of extracting what we need of course
                match = sub_pattern.match(sub_id)

                if match:
                    # Get submitter three letter code
                    name_code = match.group(1)
                    # Get institute three letter code
                    institute_code = match.group(2)
                    # Get random component of submission ID
                    random_id = match.group(3)

                    # Generate unique key, which is just $sub_id minus the final -INT
                    key = f'{name_code}-{institute_code}-{random_id}'

                    # The value we store here is the value we know we need for the SQL statement below
                    unique_submissions[key] = (key, name_code, institute_code, random_id)
                else:
                    # Warn if expected pattern doesn't match
                    print(f'Warning f{sub_id} doesn\'t conform to pattern')

            # Execute SQL statement to insert all unique submissions
            conn.executemany('INSERT INTO submissions (submission_id, name_code, institute_code, random_id) VALUES(?,?,?,?)', unique_submissions.values())

    def populate_compounds_table(self, all_data_file: Path):
        """Populate compounds table using $all_data_file
        """

        # Get database connection
        conn = self.get_conn()

        # Compile submission ID pattern
        sub_pattern = re.compile('(\w{3})-(\w{3})-(\w{8})-(\d+)')

        # Initialise compound list
        compounds = []
        compounds_noIC50 = []

        # Open the data file
        with open(all_data_file, 'r', encoding="utf8") as fh:
            # Initialise a CSV reader
            reader = csv.DictReader(fh, delimiter=',')

            # Used to identify duplicate CIDS
            cmp_idx = {}

            # Iterate rows in data file
            for cols in reader:
                # Extract fields we require
                smiles = cols['SMILES']
                comp_id = cols['CID']
                made = cols['MADE']
                inchi_key = cols['InChIKey']
                mw = cols['MW']
                r_avg_IC50 = cols['r_avg_IC50']

                # Identify and skup duplicate CIDs
                if comp_id in cmp_idx:
                    print(f'Compound ID {comp_id} already uploaded')
                    continue

                # Remember that we have seen this CID
                cmp_idx[comp_id] = 1

                # Apply pattern to submission ID
                match = sub_pattern.match(comp_id)

                if match:
                    name_code = match.group(1)
                    institute_code = match.group(2)
                    random_id = match.group(3)

                    sub_id = f'{name_code}-{institute_code}-{random_id}'

                    # Store all the data we need to upload this compound
                    if r_avg_IC50:
                        compounds.append(
                            (comp_id, sub_id, smiles, made, inchi_key, mw, r_avg_IC50))
                    else:
                        compounds_noIC50.append(
                            (comp_id, sub_id, smiles, made, inchi_key, mw))

            # Upload compounds
            conn.executemany(
                'INSERT INTO compounds (compound_id, submission_id, smiles, made, inchi_key, mw, r_avg_IC50) VALUES(?,?,?,?,?,?,?)', compounds)
            conn.executemany(
                'INSERT INTO compounds (compound_id, submission_id, smiles, made, inchi_key, mw) VALUES(?,?,?,?,?,?)', compounds_noIC50)

    def populate_assays_table(self, all_data_file: Path):
        """Populate the table submissions by reading out all of the unique submissions from $all_data_file
        """

        # Get database connection
        conn = self.get_conn()

        # Compile REGEX for submission ID
        sub_pattern = re.compile('(\w{3})-(\w{3})-(\w{8})-(\d+)')

        # Open the data file
        with open(all_data_file, 'r', encoding="utf8") as fh:
            # Initialise a CSV reader
            reader = csv.DictReader(fh, delimiter=',')

            # Create dictionary to generate unique list with
            unique_assays = {}

            # Iterate rows in data file
            for cols in reader:
                # Get all attributes
                sub_id = cols["CID"]
                f_avg_IC50 = cols["f_avg_IC50"]
                r_avg_IC50 = cols["r_avg_IC50"]
                covalent_fragment = cols["Covalent Fragment"]
                covalent_warhead = cols["covalent_warhead"]

                # The value we store here is the value we know we need for the SQL statement below
                unique_assays[sub_id] = (sub_id, r_avg_IC50, f_avg_IC50, covalent_fragment, covalent_warhead)

            # Execute SQL statement to insert all unique submissions
            conn.executemany('INSERT INTO assays (compound_id, r_avg_IC50, f_avg_IC50, covalent_fragment, covalent_warhead) VALUES(?,?,?,?,?)', unique_assays.values())

    def populate_all_tables(self, file: Path):
        self.populate_compounds_table(file)
        self.populate_assays_table(file)

    def select_compounds(self, n: int):
        """Select the SMILES string and assay data for a certain number of compounds
        """
        conn = self.get_conn()

        # Execute SQL
        data = conn.execute(f'''
            SELECT
                a.compound_id,
                a.smiles,
                b.r_avg_IC50,
                b.f_avg_IC50
            FROM
                compounds a,
                assays b
            WHERE
                a.compound_id=b.compound_id AND
                b.r_avg_IC50 <> '' AND
                b.r_avg_IC50 < 99 AND
                b.f_avg_IC50 <> '' AND
                b.f_avg_IC50 < 99
            LIMIT {n}
        ''')
        output_dict = {}
        # Print out results
        for row in data:
            print(f"Compound ID => {row[0]}, Mass Spec IC50 => {row[2]:.4f}uM, "
                  f"Fluorescence IC50 => {row[3]:.4f}uM, SMILES => {row[1]}")
            output_dict[row[0]] = row[1:]
        return output_dict


if __name__ == '__main__':
    file = Path('covid_submissions_all_info.csv')

    manager = DatabaseManager(database_path='moonshot.db')
    manager.drop_all()
    manager.create()
    manager.populate_all_tables(file)
    manager.select_compounds(10)
