install :
	gcc -g -Wall -fPIC --shared -o fp64.dll ../lua-5.3.5/src/lua53.dll  -I../lua-5.3.5/src