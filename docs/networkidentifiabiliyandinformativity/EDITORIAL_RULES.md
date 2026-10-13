# Editorial Rules

This file is the shared, user-editable source of exposition, terminology, and notation requirements for the research paper. Add new rules or amend existing rules here. These rules apply to all research prose in this workspace unless a prompt explicitly overrides them.

## Audience and exposition

1. Write for a reader with the mathematical background normally covered in a master's degree in engineering.
2. Assume elementary finite sets, conditional probability, expectations, vectors and matrices, derivatives, and constrained optimization.
3. Define specialized concepts from decision theory, game theory, partial observation, reinforcement learning, imperfect-information search, and related fields before their first use.
4. Introduce concepts in dependency order. Do not use a concept to define another concept before the first concept has been introduced.
5. Explain the purpose of a construction and what is optimized, by whom, and over which variables whenever those points are not immediate.
6. Explain important formulas and the roles of their terms without defining mathematical background already assumed above.
7. Use ordinary language and the simplest formulation that is mathematically accurate.
8. Do not introduce concepts prematurely merely because they will be useful later. Introduce them where they first perform a necessary role.
9. Avoid editorial, defensive, or meta-level commentary in the paper, including explanations of why unused notation was omitted or statements that a formula is not a definition.

## Definitions and terminology

1. Define every specialized term before its first occurrence, including broad organizing terms such as `system`, `process`, `model`, or `environment` when they are used technically rather than in their ordinary-language sense.
2. Before replacing a term, verify that the replacement is already defined or define it at the point where it is first needed.
3. Reuse an existing concept when the intended object is the same. Introduce a new term only for a genuinely different object.
4. State explicitly how each specialization, extension, or alternative representation relates to definitions introduced earlier.
5. Explain distinctions that affect the mathematics or interpretation. Eliminate distinctions that do not serve a necessary purpose.
6. Do not silently change the meaning of a term between sections.

## Symbols and notation

1. Define every symbol, function, set, index, and role before its first use.
2. Reuse established notation for the same object.
3. Introduce new notation only for a genuinely different object.
4. Never silently redefine or overload a symbol.
5. When notation is changed, audit every occurrence in the affected section and all later dependent sections.
6. State where every probability, utility, payoff, reward, score, value, and objective comes from.

## Probability notation

1. Use `\mathbb{P}` for probabilities of events and conditional probabilities of events. Do not use `\Pr`.
2. Use lowercase symbols for probability mass functions, densities, or kernels, and give distinct functions distinct names according to their roles.
3. Reserve `p` for the state-transition probability function or kernel.
4. Do not introduce an auxiliary probability function when the probability can be written clearly from already defined distributions and policies.
5. Distinguish a probability measure or event probability, written with `\mathbb{P}`, from a probability mass function or kernel, written with an appropriate lowercase symbol.

## Editing procedure

1. Read the full affected subsection before editing it.
2. Trace every new or changed term to its first occurrence and definition.
3. Trace every new or changed symbol through its later uses.
4. Do not make an isolated synonym substitution without checking whether the replacement satisfies the definition and consistency rules.
5. If the user says `do not edit`, answer the question without modifying files.
6. Keep edits focused on the requested issue unless a broader correction is required to preserve correctness or consistency.
7. Do not change established notation, terminology, structure, or style merely to make it more uniform, more elaborate, or globally unique when the existing usage is clear in context.
8. Prefer the simplest conventional notation. Do not replace a simple symbol by decorated notation such as a calligraphic symbol unless the distinction is mathematically necessary.
9. Do not make a broader change that the user did not request without first explaining why it is necessary and asking for confirmation, except when the broader change is required to keep the requested edit mathematically or syntactically correct.

## Completion audit

Before declaring an editing task complete:

1. Read the affected material linearly as a newcomer.
2. Check that every specialized term is defined before first use.
3. Check that every symbol, function, set, index, and role is defined before first use.
4. Check dependency order and explicit connections to earlier definitions.
5. Check global consistency of changed terminology and notation.
6. Check the probability-notation rules, including the absence of `\Pr` and the reserved role of `p`.
7. Check that objectives state what is optimized, by whom, and over what variables where needed.
8. Compile the LaTeX document when files were edited.
9. Review compilation warnings and report unresolved warnings or inconsistencies rather than claiming clean completion.