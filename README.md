# bbsStanBayes
early attempt to translate bbsBayes models to Stan

Initially working through the existing BBS models and providing Stan versions
Not worrying too much about keeping the specific priors consistent. Hierarchical structures will be the same, but since moving from BUGS language to Stan many of the precision priors can be better expressed as priors on SD values. Also, the original bbs models include some rather silly "uninformative" priors that place significant prior density in regions of the prameter space that are extremely unlikely.

