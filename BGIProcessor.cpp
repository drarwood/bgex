/*
   BGEN Extraction Tool (bgex)
   BGIProcessor.cpp
   Functions related to finding byte addresses of variants in BGEN files
   Also maintains user defined chr:pos:a1:a2
   Written by Andrew R Wood
   a.r.wood@exeter.ac.uk
*/

#include "BGIProcessor.h"
#include "sqlite/sqlite3.h"

using namespace std;

map<uint64_t,string> var_start_bytes;

static int callback(void *data, int argc, char **argv, char **azColName){

    // This function is called when return each row from a select query
    // typedef int (*sqlite3_callback)(
    //   void*,    // Data provided in the 4th argument of sqlite3_exec() */
    //   int,      // The number of columns in row */
    //   char**,   // An array of strings representing fields in the row */
    //   char**    // An array of strings representing column names */
    // );

    string chr = argv[0];
    string pos = argv[1];
    string a1  = argv[4];
    string a2  = argv[5];

    // remove annoying leading zeros from chr identifier
    if (chr.substr(0,1) == "0") {
        chr = chr.substr(1,1);
    }

    string id  = chr + ":" + pos  + ":" + a1 + ":" + a2;
    //cout << id << endl;
    var_start_bytes[static_cast<uint64_t>(stoul(argv[6]))] = id;
    //var_start_bytes.push_back(static_cast<uint64_t>(stoul(argv[6])));
    return 0;

}


// vector<uint64_t>& GetListofVariantStartBytes(string& f, map<string,VARIANT>& m) {
map<uint64_t,string>& GetListofVariantStartBytes(string& f, map<string,VARIANT>& m, string& c) {

    /* This function returns the start bytes of the data blocks associated with each variant in list
        if it exists...
       Parameters:
         f : ref to string holding bgen filename
         m : map of user defined chr:pos:a1:a2 -> VARIANT struct
         c : ref to string holding chromosome
    */

    var_start_bytes.clear();

    // make connection
    sqlite3* DB;
    int error = 0;
    string bgi = f+".bgi";

    error = sqlite3_open(bgi.c_str(), &DB);
    if (error) {
        std::cerr << "Error open DB " << sqlite3_errmsg(DB) << std::endl;
        exit(0);
    }


    // cycle through the variants in the list for current query to build up query string
    string where = "";
    for (map<string,VARIANT>::iterator it = m.begin(); it!=m.end(); ++it) {
        if (it->second.chr == c) {
            where = where + " (position=" + it->second.pos + " AND ((allele1='" + it->second.a1 + "' AND allele2='" + it->second.a2 + "')  OR (allele1='" + it->second.a2 + "' AND allele2='" + it->second.a1 + "'))) OR";
        }
    }
    where = where.substr(0, where.length() - 3);
    //cout << where << endl;

    // execute and call on callback
    string sql("SELECT * FROM Variant WHERE " + where +";");
    int rc = sqlite3_exec(DB, sql.c_str(), callback, 0 , NULL);

    // check if issues
    if (rc != SQLITE_OK)
        cerr << "Error SELECT" << endl;

    // close the connection
    sqlite3_close(DB);

    return var_start_bytes;

}



void MapToUserVarIDs(map<string,VARIANT>& v, map<uint64_t, string>& b) {
    /*
        This function maps bgen ids to user ids - needed if generating GRS
        Parameters:
          v: reference to map holding chr -> vector of VARIANT structs for variants on respecitve chromosome
          b: reference to map holding bgen_id -> byte start

    */
    // cycle through bgen byte starts to obtain bgen based id based on which is the dosage increasing allele

    // cycle through the map containing byte_start -> bgen variant ID has been defined in the bgi file
    for (map<uint64_t,string>::iterator it = b.begin(); it!=b.end(); ++it) {

        // assume bgi does not pre-append "chr" to id. Alter if needed to remove
        if (it->second.find("chr") != std::string::npos) {
            it->second = it->second.substr(3, it->second.length()-3);
        }

        // if there are no instances of the bgi-defined markername in the map of user defined chr:pos:a1:a2 data 
        if (!v.count(it->second)) {
            vector<string> d;
            split(it->second, ":" ,d);
            it->second = d[0]+":"+d[1]+":"+d[3]+":"+d[2];
        }
    }
}

