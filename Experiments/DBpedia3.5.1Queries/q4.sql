PREFIX  db:   <http://dbpedia.org/ontology/>
PREFIX  foaf: <http://xmlns.com/foaf/0.1/>
PREFIX  property: <http://dbpedia.org/property/>

SELECT  *
WHERE
  { 
    ?starring <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://dbpedia.org/ontology/Film>.
    ?starring db:starring ?musician .
    ?musician db:activeYearsStartYear ?activeyearsstartyear .
    ?musician db:associatedBand ?associatedband .
    ?associatedband <http://xmlns.com/foaf/0.1/page> ?page.
    ?associatedband <http://dbpedia.org/ontology/activeYearsStartYear> ?activeYearsStartYear2.
  }order by desc(?activeyearsstartyear+?activeYearsStartYear2)
LIMIT 5