PREFIX  dc:   <http://purl.org/dc/elements/1.1/>
PREFIX  :     <http://dbpedia.org/resource/>
PREFIX  rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX  dbpedia2: <http://dbpedia.org/property/>
PREFIX  geo:  <http://www.w3.org/2003/01/geo/wgs84_pos#>
PREFIX  foaf: <http://xmlns.com/foaf/0.1/>
PREFIX  owl:  <http://www.w3.org/2002/07/owl#>
PREFIX  xsd:  <http://www.w3.org/2001/XMLSchema#>
PREFIX  dbpedia: <http://dbpedia.org/>
PREFIX  rdf:  <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
PREFIX  skos: <http://www.w3.org/2004/02/skos/core#>

SELECT  ?subject ?lat ?long ?population 
WHERE
  { 
    ?area <http://www.w3.org/2004/02/skos/core#broader>  <http://dbpedia.org/resource/Category:Towns_in_England_by_county> .
    ?subject skos:subject ?area .
    ?subject geo:lat ?lat .
    ?subject geo:long ?long .
    ?subject dbpedia2:population ?population
  } order by (?lat+?long+ 0.001*?population)
LIMIT  5