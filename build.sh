mkdir -p $PREFIX/MetaLAFFA
cp -r bin $PREFIX
cp -r config $PREFIX/lib/python*/
cp lib/* $PREFIX/lib/python*/
cp -r config $PREFIX/MetaLAFFA
cp lib/* $PREFIX/MetaLAFFA
cp -r src $PREFIX/MetaLAFFA
cp Snakefile $PREFIX/MetaLAFFA
cp pipeline_steps.txt $PREFIX/MetaLAFFA
cp MetaLAFFA.py $PREFIX/MetaLAFFA
sed -i "s|{INSTALL_DIR}|$PREFIX/MetaLAFFA|" $PREFIX/lib/python*/config/file_organization.py
sed -i "s|{INSTALL_DIR}|$PREFIX/MetaLAFFA|" $PREFIX/MetaLAFFA/config/file_organization.py
sed -i "s|{CONDA_ENV}|$CONDA_DEFAULT_ENV|" $PREFIX/lib/python*/config/operation.py
sed -i "s|{CONDA_ENV}|$CONDA_DEFAULT_ENV|" $PREFIX/MetaLAFFA/config/operation.py