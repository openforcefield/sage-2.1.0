# remove parameterize keyword from all lines in FF
# Author: Pavan Behara
import click


@click.command()
@click.option(
    "-ff", "--ff", "ff_to_modify", type=click.STRING, default="force-field.offxml"
)
def main(ff_to_modify):
    from openff.toolkit.typing.engines.smirnoff import ForceField

    forcefield = ForceField(ff_to_modify, allow_cosmetic_attributes=True)
    for handler in forcefield.registered_parameter_handlers:
        print(handler)
        if len(forcefield[handler].parameters) == 0:
            continue
        print(forcefield[handler].parameters)
        for par in forcefield[handler].parameters:
            if hasattr(par, "_parameterize"):
                par.delete_cosmetic_attribute("parameterize")

    forcefield.to_file("force-field-no-cosmetic-params.offxml")


if __name__ == "__main__":
    main()
