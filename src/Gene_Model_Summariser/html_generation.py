#pillar 3 - HTML Report Generation

from pathlib import Path
from jinja2 import Environment, FileSystemLoader, select_autoescape


def render_report(context: dict) -> str:
    # 1) Find the templates folder (relative to this Python file)
    templates_folder = Path(__file__).resolve().parent / "templates"

    # 2) Create a Jinja environment that knows where templates live
    env = Environment(loader=FileSystemLoader(str(templates_folder)), autoescape=select_autoescape(["html", "xml"]),)
    # 3) Load your template file by name
    template = env.get_template("report_template.html")

    # 4) Render the template using your context dictionary
    html_output = template.render(**context)

    return html_output
