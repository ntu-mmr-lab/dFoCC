# Cartesian product -- to reduce the use nested loops
from itertools import product

# for type hint
from typing import Callable, NamedTuple, cast

from Bio.PDB.Structure import Structure

ParamTuple = tuple[int, ...]
RCModFunction = Callable[
    [Structure, ParamTuple],
    None
]
PropFunction = Callable[
    [Structure, ParamTuple],
    tuple[bool, dict[str, float]]
]

class ReactionCoordinate:
    def __init__(
        self,
        rc_mod_fn: RCModFunction,
        init_step_size: ParamTuple,
        step_per_run: ParamTuple | None = None,
        min_step: ParamTuple | None = None,
    ):
        if not step_per_run:
            step_per_run = cast(ParamTuple, (5, ) * len(init_step_size))
        if not min_step:
            min_step = cast(ParamTuple, (1, ) * len(init_step_size))
        assert len(init_step_size) == len(step_per_run) == len(min_step), \
            f'The Reaction Coordinate does not have a consistent set of parameter: {len(init_step_size) = }, {len(step_per_run) = }, {len(min_step) = }'
        self.rc_mod_fn = rc_mod_fn
        self.init_step_size = init_step_size
        self.step_per_run = step_per_run
        self.min_step = min_step
        self.accept_params = len(init_step_size)

    def __call__(self, struct: Structure, params: ParamTuple):
        assert len(params) == self.accept_params, \
            f'The paramter size given ({len(params)}) is not the same as the size ({self.accept_params}) required for this reaction coordinate ({self.__class__.__name__})'
        return self.rc_mod_fn(struct, params)

class ReactionCoordinateSeries:
    class ParametersReturnType(NamedTuple):
        is_last_round: bool
        params_list: list[tuple[int, ...]]

    def __init__(self):
        self.rc_list: list[ReactionCoordinate] = []
        self.property_fn_list: list[PropFunction] = []

    def append_rc(self, rc: ReactionCoordinate):
        self.rc_list.append(rc)

    def append_properties(self, prop_fn: PropFunction):
        self.property_fn_list.append(prop_fn)

    def use_rc(
        self,
        init_step_size: ParamTuple,
        step_per_run: ParamTuple | None = None,
        min_step: ParamTuple | None = None,
        rc_mod_fn: RCModFunction | None = None
    ):
        def decorator(rc_mod_fn: RCModFunction):
            rc: ReactionCoordinate = ReactionCoordinate(rc_mod_fn, init_step_size, step_per_run, min_step)
            self.rc_list.append(rc)
        if callable(rc_mod_fn):
            decorator(rc_mod_fn)
            return decorator
        return decorator

    def get_best_params(self):
        return (0, ) * sum(rc.accept_params for rc in self.rc_list)

    def generate_params_list(self, round_num: int):
        params_it_list = []
        is_last_round = True
        for rc in self.rc_list:
            for i in range(rc.accept_params):
                # in each round, every parameter gets halved, until it reaches the minimal step
                round_step_size = max(rc.min_step[i], rc.init_step_size[i] // 2 ** (round_num - 1))
                # if odd number of steps is set, the interval is set evenly on both sides of 0
                # if even number of steps is set, the interval is bias on the positive side
                # no matter how, the iterator must include 0
                range_start = -round_step_size * ((rc.step_per_run[i] - 1) // 2)
                range_end = round_step_size * (rc.step_per_run[i] // 2 + 1)
                params_it_list.append(range(range_start, range_end, round_step_size))
                # if any one of the parameter does not reach the minimal step, it is not the last round
                if round_step_size != rc.min_step[i]:
                    is_last_round = False

        # confirm the best structure (all 0's) is included in the loop
        # otherwise, the program below cannot compare the results with the last cycle
        params_list = list(product(*params_it_list))
        best_params = self.get_best_params()
        assert best_params in params_list, f'The best parameters in the last round / the default parameters {best_params} is not included in the parameter list'

        return self.ParametersReturnType(
            is_last_round,
            params_list # lambda: product(params_it)
        )

    def modify_moiety(self, structure: Structure, params: tuple[int, ...]):
        assert len(params) == sum(rc.accept_params for rc in self.rc_list), \
            f'The parameter size ({len(params)}) given is not the same as the size required for all reaction coordinates ({" + ".join(str(rc.accept_params) for rc in self.rc_list)})'
        params_count = 0
        for rc in self.rc_list:
            rc(structure, params[params_count:params_count + rc.accept_params])
            params_count += rc.accept_params

    def calculate_moiety_properties(self, structure: Structure, params: ParamTuple):
        is_ok = True
        property_dict: dict[str, float] = {}
        for property_fn in self.property_fn_list:
            is_this_ok, this_property_dict = property_fn(structure, params)
            if not is_this_ok:
                is_ok = False
            property_dict.update(this_property_dict)
        return is_ok, property_dict

# if __name__ == '__main__':
#     rcs = ReactionCoordinateSeries()
#
#     @rcs.use_rc(init_step_size=(3, 30))
#     def first_rc(_: Structure, params):
#         print(params)
#
#     @rcs.use_rc(init_step_size=(20, ))
#     def second_rc(_: Structure, params):
#         print(params)
#
#     rcs.use_rc(
#         rc_mod_fn=lambda _, params: print(params),
#         init_step_size=(3, ),
#         step_per_run=(1, ),
#     )
#
#     rcs.append_rc(ReactionCoordinate(
#         rc_mod_fn=lambda _, params: print(params),
#         init_step_size=(2, 10),
#         step_per_run=(1, 1),
#     ))
#
#     print(rcs.generate_params_list(3))
#     print(rcs.generate_params_list(100))
#     rcs.modify_moiety(None, (3, 3, 3, 10, 12, 14)) # type: ignore
#
#     @rcs.append_properties
#     def calc_props(_):
#         return False, { 'test': 0. }
#
#     rcs.append_properties(lambda _: (True, { 'str': 1.1 }))
#
#     print(rcs.calculate_moiety_properties(None)) # type: ignore
#
