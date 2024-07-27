PREFIX  resource: <http://dbpedia.org/resource/>
PREFIX  property: <http://dbpedia.org/property/>

SELECT DISTINCT  ?id ?officialName ?populationTotal
WHERE
  { ?id property:settlementType resource:City .
    ?id property:officialName ?officialName .
    ?id property:populationTotal ?populationTotal.
    } order by ?populationTotal limit 5