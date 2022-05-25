# `UNM_BIOCOMP_DEPICT`

Molecular depiction classes, including mol2img servlets and depict webapp.
Also includes 2D alignment tools (formerly in `unm_biocomp_molalign`).

## Dependencies

* Java 8
* Maven 3.5+
* ChemAxon (19.3.0)
* Access to [ChemAxon Maven repository](https://hub.chemaxon.com)
(see [documentation](https://docs.chemaxon.com/display/docs/Public+Repository))
  * Requires API key.
* `unm_biocomp_util`

## Compiling

```
mvn clean install
```

## Usage

```
mvn exec:java -Dexec.mainClass="edu.unm.health.biocomp.depict.mol2img_app" -Dexec.args="-i quinine.smi -o quinine.png"
```

