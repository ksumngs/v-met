import nf_core.schema

def param_parse(app, config):
    schema = nf_core.schema.PipelineSchema()
    schema.get_schema_path("nextflow_schema.json")
    schema.load_schema()
    docs = schema.print_documentation(
        None,
        "markdown",
        False,
        ["parameter", "description", "type", "default", "required"]
    ).splitlines()
    docs[0] = "# Complete Parameter Reference"
    docs[2] = "Every single parameter in the entire pipeline is described here. If you need to tweak every little detail of how your samples are analyzed, you've come to the right place! Note that *every* parameter is described here, even those that shouldn't be used, so proceed with caution!"

    with open("docs/parameters.md", "w") as f:
        for l in docs:
            f.write(f"{l}\n")

def setup(app):
    app.connect('config-inited', param_parse)
