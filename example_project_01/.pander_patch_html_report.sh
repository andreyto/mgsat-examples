#!/bin/bash
find . -name '*.html' | xargs sed -i'.sed.bak' -e 's/href="http:/href="https:/'
find . -name '*.html' | xargs sed -i'.sed.bak' -e 's/src="http:/src="https:/'
find . -name '*.html' | xargs sed -i'.sed.bak' -e '/stylesheets\/skeleton.css/s/$/\'$'
''<style>\'$'
''\.container \{ width: 100\%; \}\'$'
''\.container \.twelve\.columns \{   width: 80\%; \}\'$'
''\.container \.three\.columns \{   width: 20\%; \}\'$'
''<\/style>\'$'
''/'
rm *.sed.bak

