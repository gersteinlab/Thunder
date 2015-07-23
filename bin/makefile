all: 
	cd ~/Programming/JavaProjects/Thunder/bin
	jar cmf META-INF/MANIFEST.MF Thunder.jar *
	jar uf Thunder.jar -C ~/Programming/JAR-files/commons-io-2.4/ org
	jar uf Thunder.jar -C ~/Programming/JAR-files/commons-cli-1.2/ org
	jar uf Thunder.jar -C ~/Programming/JAR-files/commons-lang3-3.3.2/ org
	jar uf Thunder.jar -C ~/Programming/JAR-files/commons-math3-3.5/ org
	jar uf Thunder.jar -C ~/Programming/JAR-files/sqlite-jdbc-3.7.2/ org
	jar uf Thunder.jar -C ~/Programming/JAR-files/sqlite-jdbc-3.7.2/ native
	jar uf Thunder.jar -C ~/Programming/JAR-files/sam-1.96/ net
	java -jar Thunder.jar
